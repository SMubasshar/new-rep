{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "import time\n",
    "\n",
    "# Install FEniCS\n",
    "try:\n",
    "    import dolfin\n",
    "except ImportError as e:\n",
    "    !apt-get install -y -qq software-properties-common\n",
    "    !add-apt-repository -y ppa:fenics-packages/fenics\n",
    "    !apt-get update -qq\n",
    "    !apt install -y --no-install-recommends fenics\n",
    "    !sed -i \"s|#if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 8 && PETSC_VERSION_RELEASE == 1|#if 1|\" /usr/include/dolfin/la/PETScLUSolver.h\n",
    "    !rm -rf /usr/lib/python3/dist-packages/mpi4py*\n",
    "    !rm -rf /usr/lib/python3/dist-packages/petsc4py*\n",
    "    !rm -rf /usr/lib/python3/dist-packages/slepc4py*\n",
    "    !rm -rf /usr/lib/petsc/lib/python3/dist-packages/dolfin*\n",
    "    !rm -rf /usr/lib/petsc/lib/python3/dist-packages/mshr*\n",
    "    !wget \"https://drive.google.com/uc?export=download&id=1cT_QBJCOW_eL3BThnval3bcpb8o0w-Ad\" -O /tmp/mpi4py-2.0.0-cp37-cp37m-linux_x86_64.whl\n",
    "    !wget \"https://drive.google.com/uc?export=download&id=119i49bxlGn1mrnhTNmOvM4BqmjrT9Ppr\" -O /tmp/petsc4py-3.7.0-cp37-cp37m-linux_x86_64.whl\n",
    "    !wget \"https://drive.google.com/uc?export=download&id=1-1tVfu8qz3bRC2zvR8n3RESpesWqNnn6\" -O /tmp/slepc4py-3.7.0-cp37-cp37m-linux_x86_64.whl\n",
    "    !wget \"https://drive.google.com/uc?export=download&id=1-3qY4VIJQaXVO1HfGQIzTIURIeJbvX-9\" -O /tmp/fenics_dolfin-2019.2.0.dev0-cp37-cp37m-linux_x86_64.whl\n",
    "    !wget \"https://drive.google.com/uc?export=download&id=1-5SMjgjMuee_9WLeYtGe8N_lvipWEN7W\" -O /tmp/mshr-2019.2.0.dev0-cp37-cp37m-linux_x86_64.whl\n",
    "    !pip3 install /tmp/mpi4py-2.0.0-cp37-cp37m-linux_x86_64.whl --upgrade\n",
    "    !pip3 install /tmp/petsc4py-3.7.0-cp37-cp37m-linux_x86_64.whl --upgrade\n",
    "    !pip3 install /tmp/slepc4py-3.7.0-cp37-cp37m-linux_x86_64.whl --upgrade\n",
    "    !pip3 install /tmp/fenics_dolfin-2019.2.0.dev0-cp37-cp37m-linux_x86_64.whl --upgrade\n",
    "    !pip3 install /tmp/mshr-2019.2.0.dev0-cp37-cp37m-linux_x86_64.whl --upgrade\n",
    "    !pip3 -q install --upgrade sympy\n",
    "    import dolfin\n",
    "\n",
    "from dolfin import *; from mshr import *\n",
    "\n",
    "import dolfin.common.plotting as fenicsplot\n",
    "\n",
    "from matplotlib import pyplot as plt"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.10"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
