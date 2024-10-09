import numpy as np
from lsst.ts.wep.estimation import WfAlgorithmFactory, WfEstimator
from lsst.ts.wep.image import Image
from lsst.ts.wep.task.estimateZernikesBase import EstimateZernikesBaseConfig
from lsst.ts.wep.task.estimateZernikesDanishTask import EstimateZernikesDanishConfig
from lsst.ts.wep.task.estimateZernikesTieTask import EstimateZernikesTieConfig
from lsst.ts.wep.utils import (
    WfAlgorithmName,
    getCameraFromButlerName,
    getTaskInstrument,
)

from lsst.daf import butler as dafButler


def stack_donut_wep_im_refactor(
    donut_stamps,
    n=2,
    indices=None,
    pixel_stack="mean",
    use_mask=False,
    replace_zeros_with_nans=True,
    after_avg_fill_with_bkgnd=True,
    use_mean_bkgnd=True,
    use_range_bkgnd=False,
    bkgnd_factor=10,
):

    stacked = {}

    if indices is None:
        indices = np.arange(n)

    else:
        n = len(indices)

    stacked["N"] = n
    stacked["indices"] = indices

    # create updated donut stamps consisting of the stacked images

    # initialize empty arrays of the correct dimension
    image_stack = []
    bkgnd_stack = []

    # add donut image arrays, store the x,y
    xs = []
    ys = []

    j = 0
    for i in indices:
        # print(f'stacking {i}')
        stamp = donut_stamps[i]
        image = stamp.stamp_im.image.array

        if use_mask:
            # has values like 0, 1024,  2048 etc
            mask = stamp.stamp_im.mask.array
            binary_mask = np.array(mask > 0).astype(int)  # has values 0 or 1
            arr = image * binary_mask
            if replace_zeros_with_nans:
                arr[arr == 0] = np.nan  # replace zeros with NaNs
        else:
            arr = image

        image_stack.append(arr)
        bkgnd_stack.append(image)

        # that already transposes each field angle ...
        fieldxy = stamp.calcFieldXY()[::-1]
        xs.append(fieldxy[0])
        ys.append(fieldxy[1])

        j += 1

        # print('stacked: ', indices)
        # store the original fieldXY
        stacked["fieldX"] = xs
        stacked["fieldY"] = ys

    # store the stacked image
    if pixel_stack == "sum":
        donut_stacked_array = np.sum(image_stack, axis=0)
    elif pixel_stack == "mean":
        donut_stacked_array = np.mean(image_stack, axis=0)
    elif pixel_stack == "nanmean":
        donut_stacked_array = np.nanmean(image_stack, axis=0)
        bkgnd_mean = np.nanmean(bkgnd_stack, axis=0)

    # replace np.nan with mean background...
    if use_mask and after_avg_fill_with_bkgnd:
        mask_nans = np.isnan(donut_stacked_array)
        if use_mean_bkgnd:
            donut_stacked_array[mask_nans] = bkgnd_mean[mask_nans]
            stacked["bkgnd"] = bkgnd_mean
        elif use_range_bkgnd:

            # change interval (0,1) to (-0.5,0.5) , and increase by a
            # factor to very small value or larger one
            size = np.shape(donut_stacked_array)
            bkgnd = bkgnd_factor * (np.random.random_sample(size=size) - 0.5)
            print(
                f"adding background in range {min(np.ravel(bkgnd))}\
to {max(np.ravel(bkgnd))}"
            )
            stacked["bkgnd"] = bkgnd
            donut_stacked_array[mask_nans] = bkgnd[mask_nans]

    stacked["donutStackedArray"] = donut_stacked_array

    # find the midpoint fieldXY coordinates
    # average in the same way x and y coordinates
    # for the intra and extra donut alike
    # treat them as already transposed
    stacked["fieldXmean"] = np.mean(xs)
    stacked["fieldYmean"] = np.mean(ys)

    # store information on how pixels were combined
    stacked["pixelStack"] = pixel_stack

    # part of _setWepImage(self):
    """Return a ts.wep.image.Image object for the stamp.

    Note that the information from the butler is in the data visualization
    coordinate system (DVCS), but the WEP Image is in the global camera
    coordinate system (CCS). These coordinate systems are related by a
    transpose. See sitcomtn-003.lsst.io for more information.

    Furthermore, CWFS images that arrive from the butler are rotated with
    respect to the science sensors. The info in the WEP Images has been
    de-rotated so that everything aligns with the global coordinate system
    used by the science sensors.
    """
    # both the camera name and the detector would be all
    # identical for all stamps in here ,
    # hence it would be the same Euler angle

    # stamp.cam_name  yields eg. LSSTCam;
    # taking it from the stamp
    # allows to stack also images from LSSTComCam,  or FAM ...
    camera = getCameraFromButlerName(stamp.cam_name)

    detector = camera.get(stamp.detector_name)

    # Get the rotation with respect to the science sensors
    eulerz = -detector.getOrientation().getYaw().asDegrees()
    nrot = int(eulerz // 90)
    if not np.isclose(eulerz % 90, 0):
        raise RuntimeError(
            f"The detector is rotated {-eulerz} deg with respect to the science "
            "sensors, but _setWepImage() only works for sensors whose rotations "
            "are an integer multiple of 90 deg."
        )

    # Rotate to orientation of science sensors
    image = np.rot90(donut_stacked_array, nrot)

    # Transpose the image (DVCS -> CCS)
    image = image.T

    # Get the field angle, and transpose (DVCS -> CCS)
    # fieldAngle = stamp.calcFieldXY()
    # the x,y are already a transpose of calcFieldXY() above
    fieldangle = (stacked["fieldXmean"], stacked["fieldYmean"])

    # Determine the blend offsets
    if stamp.blend_centroid_positions.size > 0:
        # Get the offsets in the original pixel coordinates
        blendoffsets = stamp.blend_centroid_positions - stamp.centroid_position

        # Rotate the coordinates (by -90 each time)
        # to match the science sensors
        rotmat = np.array([[0, 1], [-1, 0]])
        if stamp.defocal_type == "extra":
            rotmat = np.linalg.matrix_power(rotmat, nrot + 2)
        else:
            rotmat = np.linalg.matrix_power(rotmat, nrot)
        blendoffsets = np.transpose(rotmat @ blendoffsets.T)

    else:
        blendoffsets = None

    # Package everything in an Image object
    wepimage = Image(
        image=image,
        fieldAngle=fieldangle,
        defocalType=stamp.defocal_type,
        bandLabel=stamp.bandpass,
        blendOffsets=blendoffsets,
    )

    stacked["wep_im"] = wepimage

    return stacked


# get the data directly from the butler
butler_root_path = "/sdf/group/rubin/shared/scichris/DM-42718_WET-006/newRepo"
output_collection = "all_states_direct_stamps"

butler = dafButler.Butler(butler_root_path)
registry = butler.registry

dataset_refs = registry.queryDatasets(
    "donutStampsExtra",
    collections=[output_collection],
    where=f"instrument='LSSTCam' and detector.purpose='WAVEFRONT' ",
).expanded()

print(len(list(dataset_refs)))
refs = []
for ref in dataset_refs:
    refs.append(ref)


for use_mask in [True, False]:
    for method in ["tie", "danish"]:
        add_bkgnd = False

        # set the file name title
        string = "no"
        if use_mask:
            string = "use"

        bkgnd = "no"
        if add_bkgnd:
            bkgnd = "with"
        fname = f"wep_direct_stacking_{string}_mask_{method}_{bkgnd}_bkgnd.npy"

        results = {}
        for state in range(1, 101):
            results[state] = {}

        # store per state ,  per corner ,  the intra/extra  stamps
        for ref in refs:  # [:1]:
            state = int(str(ref.dataId.visit.id)[-3:])
            raft = ref.dataId.detector.raft

            print(f"Reading state {state},  {raft} ")

            stamps_extra = butler.get(
                "donutStampsExtra", dataId=ref.dataId, collections=[output_collection]
            )
            stamps_intra = butler.get(
                "donutStampsIntra", dataId=ref.dataId, collections=[output_collection]
            )

            # Initialize WEP following calcZernikes =-> estimateZernikes

            # read the camera and detector from the extra-focal donut, since
            # the instrument would be the same for all
            cam_name = stamps_extra[0].cam_name
            detector_name = stamps_extra[0].detector_name

            instrument = getTaskInstrument(
                cam_name,
                detector_name,
                None,
            )

            estimateZernikesBaseConfig = EstimateZernikesBaseConfig()

            if method == "danish":
                wfAlgoName = WfAlgorithmName.Danish
                estimateZkConfig = EstimateZernikesDanishConfig()
            elif method == "tie":
                wfAlgoName = WfAlgorithmName.TIE
                estimateZkConfig = EstimateZernikesTieConfig()

            algoConfig = {
                key: val
                for key, val in estimateZkConfig.toDict().items()
                if key not in estimateZernikesBaseConfig._fields.keys()
            }
            wfAlgoConfig = WfAlgorithmFactory.createWfAlgorithm(wfAlgoName, algoConfig)

            # Create the wavefront estimator
            wfEst = WfEstimator(
                algoName=wfAlgoName,
                algoConfig=wfAlgoConfig,
                instConfig=instrument,
                jmax=estimateZernikesBaseConfig.maxNollIndex,
                startWithIntrinsic=estimateZernikesBaseConfig.startWithIntrinsic,
                returnWfDev=estimateZernikesBaseConfig.returnWfDev,
                return4Up=estimateZernikesBaseConfig.return4Up,
                units=estimateZernikesBaseConfig.units,
                saveHistory=estimateZernikesBaseConfig.saveHistory,
            )

            # take the mean defocal offset from the first donut in each
            # detector; this would be the same for all donuts
            donut_extra = stamps_extra[0]
            donut_intra = stamps_intra[0]
            # the mean defocal offset would be the same for all considered donuts
            defocalOffset = np.mean(
                [
                    donut_extra.defocal_distance,
                    donut_intra.defocal_distance,
                ]
            )
            # print(state, raft, defocalOffset)
            wfEst.instrument.defocalOffset = defocalOffset / 1e3  # m -> mm

            # stack donuts and estimate as well
            stackedExtra = stack_donut_wep_im_refactor(
                stamps_extra,
                N=len(stamps_extra),
                pixel_stack="nanmean",
                use_mask=use_mask,
                after_avg_fill_with_bkgnd=add_bkgnd,
            )

            stackedIntra = stack_donut_wep_im_refactor(
                stamps_intra,
                N=len(stamps_intra),
                pixel_stack="nanmean",
                use_mask=use_mask,
                after_avg_fill_with_bkgnd=add_bkgnd,
            )

            zkst = wfEst.estimateZk(stackedExtra["wep_im"], stackedIntra["wep_im"])
            results[state][raft] = zkst

        np.save(fname, results, allow_pickle=True)
        print(f"Saved as {fname}")
