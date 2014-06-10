from distutils.core import setup

setup(name="ctmixtures",
      version="1.0",
      packages=['ctmixtures',
                'ctmixtures.utils',
                'ctmixtures.analysis',
                'ctmixtures.data',
                'ctmixtures.traits',
                'ctmixtures.population',
                'ctmixtures.dynamics',
                'ctmixtures.rules'],
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
