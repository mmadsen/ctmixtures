from distutils.core import setup

setup(name="ctmixtures",
      version="1.0",
      package_dir = {'': 'ctmixtures'},
      scripts = [
          'admin/ctmixtures-planner.py',
          'admin/ctmixtures-priorsampler-runbuilder.py',
          'admin/ctmixtures-export-data.py',
          'simulations/sim-ctmixture-notimeaveraging.py',
          'simulations/sim-ctmixture-timeaveraging.py'
      ],
      author='Mark E. Madsen',
      author_email='mark@madsenlab.org',
      url='https://github.com/mmadsen/ctmixtures',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: Apache Software License',
          'Operating System :: POSIX',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering',
      ]
      )
