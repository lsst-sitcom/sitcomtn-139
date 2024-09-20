# See the Documenteer docs for how to customize conf.py:
# https://documenteer.lsst.io/technotes/

from documenteer.conf.technote import *  # noqa F401 F403

extensions += ['nbsphinx']

# Assume notebook is pre-computed; we don't have the LSST Pipelines on Travis.
nbsphinx_execute = 'never'

nb_execution_mode = "off"