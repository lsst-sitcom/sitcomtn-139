.. image:: https://img.shields.io/badge/sitcomtn--139-lsst.io-brightgreen.svg
   :target: https://sitcomtn-139.lsst.io
.. image:: https://github.com/lsst-sitcom/sitcomtn-139/workflows/CI/badge.svg
   :target: https://github.com/lsst-sitcom/sitcomtn-139/actions/

##############################################
WET-006: pairing vs stacking with WET-001 data
##############################################

SITCOMTN-139
============

We compare the default approach for wavefront estimation: using individual donut pairs and averaging their final Zk estimate (pairing) to aggregating the donut stamps prior to wavefront fitting, and using their combined information to obtain a single set of Zernike polynomials per stacked donut image (stacking). 

**Links:**

- Publication URL: https://sitcomtn-139.lsst.io
- Alternative editions: https://sitcomtn-139.lsst.io/v
- GitHub repository: https://github.com/lsst-sitcom/sitcomtn-139
- Build system: https://github.com/lsst-sitcom/sitcomtn-139/actions/


Build this technical note
=========================

You can clone this repository and build the technote locally if your system has Python 3.11 or later:

.. code-block:: bash

   git clone https://github.com/lsst-sitcom/sitcomtn-139
   cd sitcomtn-139
   make init
   make html

Repeat the ``make html`` command to rebuild the technote after making changes.
If you need to delete any intermediate files for a clean build, run ``make clean``.

The built technote is located at ``_build/html/index.html``.

Publishing changes to the web
=============================

This technote is published to https://sitcomtn-139.lsst.io whenever you push changes to the ``main`` branch on GitHub.
When you push changes to a another branch, a preview of the technote is published to https://sitcomtn-139.lsst.io/v.

Editing this technical note
===========================

The main content of this technote is in ``index.rst`` (a reStructuredText file).
Metadata and configuration is in the ``technote.toml`` file.
For guidance on creating content and information about specifying metadata and configuration, see the Documenteer documentation: https://documenteer.lsst.io/technotes.
