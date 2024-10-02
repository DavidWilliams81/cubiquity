#!/usr/bin/env python3

from io import BytesIO
from urllib.request import urlopen
from zipfile import ZipFile
from tarfile import TarFile

zipurl = 'https://s3.eu-west-2.amazonaws.com/dpw81.public/cubiquity_data_v6.tar.xz'
with urlopen(zipurl) as zipresp:
    with TarFile.open(fileobj=BytesIO(zipresp.read())) as zfile:
        zfile.extractall('.')
