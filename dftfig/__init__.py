"""dftfig - automated DFT and figure pipeline for ORCA."""
import os

__version__ = "0.9.1"

# Pick a file-writing matplotlib backend before anything imports matplotlib.
#
# dftfig only ever saves figures to disk; it never opens a window. If the environment
# has PyQt in it - and chemistry environments often do, because PyMOL pulls it in -
# matplotlib may select a Qt backend on import, and then the very first plt.figure()
# blocks inside QApplication creation waiting on a display that a batch run has no use
# for. The symptom is a figures stage that prints its banner and then hangs forever at
# no CPU cost, which looks exactly like a slow render and is not one.
#
# Anyone who deliberately exports MPLBACKEND still gets their choice; an empty value
# is not a choice, and matplotlib ignores it and auto-selects, so treat it as unset.
#
# This environment variable, not the late guard in figstyle, is what actually prevents
# the hang: matplotlib consults it when it resolves the backend, which happens inside
# the first plt.figure(). By then a guard that inspects rcParams sees only an
# unresolved sentinel and cannot tell what is about to be chosen.
if not os.environ.get("MPLBACKEND"):
    os.environ["MPLBACKEND"] = "Agg"
