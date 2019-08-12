# Metadata

## Standard dcraw metadata

From a Canon CR2 file, `dcraw -i -v` will print metadata.

```
dcraw -i -v capture000003.cr2 

Filename: capture000003.cr2
Timestamp: Sat Jul 13 21:38:06 2019
Camera: Canon EOS DIGITAL REBEL XTi
Owner: unknown
ISO speed: 100
Shutter: 1/50.0 sec
Aperture: f/inf
Focal length: 0.0 mm
Embedded ICC profile: no
Number of raw images: 1
Thumb size:  3888 x 2592
Full size:   3948 x 2622
Image size:  3906 x 2602
Output size: 3906 x 2602
Raw colors: 3
Filter pattern: RG/GB
Daylight multipliers: 2.423860 0.921348 1.151113
Camera multipliers: 1997.000000 1080.000000 2333.000000 1080.000000
```

When writing TIFF format using the `-T` command line option, dcraw 
will write metada to the file. The following is the metadata written
from a Canon CR2 file, extracted using `exiv2`.

```
dcraw -4 -T capture000003.cr2 
exiv2 pr capture000003.tiff 

File name       : capture000003.tiff
File size       : 60982324 Bytes
MIME type       : image/tiff
Image size      : 3906 x 2602
Camera make     : Canon
Camera model    : EOS DIGITAL REBEL XTi
Image timestamp : 
Image number    : 
Exposure time   : 0.019999 s
Aperture        : F-2.1e+03
Exposure bias   : 
Flash           : 
Flash bias      : 
Focal length    : 0.0 mm
Subject distance: 
ISO speed       : 100
Exposure mode   : 
Metering mode   : 
Macro mode      : 
Image quality   : 
Exif Resolution : 3906 x 2602
White balance   : 
Thumbnail       : None
Copyright       : 
Exif comment    : 

```

The same TIFF file, using python3 imageio.
```
python3
Python 3.7.3 (default, Jul 12 2019, 12:37:27) 
[Clang 10.0.0 (clang-1000.11.45.5)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> import imageio
>>> im = imageio.imread('capture000003.tiff')
>>> im.meta
Dict([('is_fluoview', False), 
    ('is_nih', False), 
    ('is_micromanager', False), 
    ('is_ome', False), 
    ('is_reduced', 0), 
    ('is_sgi', False), 
    ('is_shaped', None), 
    ('is_stk', False), 
    ('is_tiled', False), 
    ('compression', <COMPRESSION.NONE: 1>), 
    ('is_mediacy', False), 
    ('description', ''), 
    ('description1', ''), 
    ('is_imagej', None), 
    ('software', 'dcraw v9.28'), 
    ('datetime', datetime.datetime(2019, 7, 13, 21, 38, 6))])
```

