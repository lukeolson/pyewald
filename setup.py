from distutils.core import setup

DESCRIPTION = 'Ewald Summation in Python'
LONG_DESCRIPTION = """pyewald: Ewald Summation in Python

Code details are found at http://github.com/lukeolson/pyewald
"""
NAME = 'pyewald'
AUTHOR = 'Luke Olson'
AUTHOR_EMAIL = 'luke.olson@gmail.com'
MAINTAINER = 'Luke Olson'
MAINTAINER_EMAIL = 'luke.olson@gmail.com'
URL = 'http://github.com/lukeolson/pyewald'
DOWNLOAD_URL = 'http://github.com/lukeolson/pyewald'
LICENSE = 'MIT'

setup(name=NAME,
      version='0.1',
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      url=URL,
      download_url=DOWNLOAD_URL,
      license=LICENSE,
      packages=['pyewald',
                'pyewald.tests'],
      )
