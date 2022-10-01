## Camera RAW Metadata

Reading and writing Camera metadata is a challenge. These notes are
more of a work in progress than a summary of reults.

### Standard dcraw metadata, CR2 format file

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

The TIFF file is relatively sane, here is an extract from EXIF metadata.
Note that the values in the dictionary are rarely anything we can simply use,
they're either in some format that should be converted or the references
to memory locations.
```
>>> for key in im.meta['EXIF_MAIN'].keys():
...     print(key, im.meta['EXIF_MAIN'][key])
... 
ExifVersion b'0221'
ComponentsConfiguration b'\x01\x02\x03\x00'
ShutterSpeedValue (369876, 65536)
DateTimeOriginal 2019:07:13 21:38:06
DateTimeDigitized 2019:07:13 21:38:06
ApertureValue (2147483648, 1)
ExposureBiasValue (0, 3)
MeteringMode 1
UserComment b'\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00'
Flash 16
FocalLength (0, 1)
ColorSpace 65535
ExifImageWidth 3888
ExifInteroperabilityOffset 5902
FocalPlaneXResolution (3888000, 877)
Make Canon
Model Canon EOS DIGITAL REBEL XTi
YCbCrCoefficients ((299, 1000), (587, 1000), (114, 1000))
Orientation 1
YCbCrPositioning 2
FocalPlaneYResolution (2592000, 582)
ExifImageHeight 2592
FocalPlaneResolutionUnit 2
XResolution (72, 1)
YResolution (72, 1)
ExposureTime (1, 50)
FNumber (0, 1)
ExposureProgram 1
CustomRendered 0
ISOSpeedRatings 100
ResolutionUnit 2
Gamma (22, 10)
ExposureMode 1
FlashPixVersion b'0100'
WhiteBalance 1
DateTime 2019:07:13 21:38:06
WhitePoint ((313, 1000), (329, 1000))
PrimaryChromaticities ((64, 100), (33, 100), (21, 100), (71, 100), (15, 100), (6, 100))
SceneCaptureType 0
ExifOffset 320

```
After some inspection of imageio's github Issues page it appears like imageio
does NOT currently provide any usable metadata handling. We may have to skip
metadata, look to another module, or fall back on system calls to a
utility like exiftool (which is not portable).

### Metadata Read By ExifRead

Using ExifRead (v2.3.2):
```
2021-03-14 15:23:06,965 | DEBUG | AstroPhotography |   key=Image ImageLength               val=1288
2021-03-14 15:23:06,965 | DEBUG | AstroPhotography |   key=Image BitsPerSample             val=[8, 8, 8]
2021-03-14 15:23:06,965 | DEBUG | AstroPhotography |   key=Image Compression               val=JPEG (old-style)
2021-03-14 15:23:06,965 | DEBUG | AstroPhotography |   key=Image Make                      val=Canon
2021-03-14 15:23:06,965 | DEBUG | AstroPhotography |   key=Image Model                     val=Canon EOS DIGITAL REBEL XTi
2021-03-14 15:23:06,965 | DEBUG | AstroPhotography |   key=Image StripOffsets              val=80168
2021-03-14 15:23:06,965 | DEBUG | AstroPhotography |   key=Image Orientation               val=Horizontal (normal)
2021-03-14 15:23:06,966 | DEBUG | AstroPhotography |   key=Image StripByteCounts           val=228169
2021-03-14 15:23:06,966 | DEBUG | AstroPhotography |   key=Image XResolution               val=72
2021-03-14 15:23:06,966 | DEBUG | AstroPhotography |   key=Image YResolution               val=72
2021-03-14 15:23:06,966 | DEBUG | AstroPhotography |   key=Image ResolutionUnit            val=Pixels/Inch
2021-03-14 15:23:06,966 | DEBUG | AstroPhotography |   key=Image DateTime                  val=2014:01:18 21:48:57
2021-03-14 15:23:06,966 | DEBUG | AstroPhotography |   key=Image ExifOffset                val=270
2021-03-14 15:23:06,966 | DEBUG | AstroPhotography |   key=Thumbnail JPEGInterchangeFormat  val=78336
2021-03-14 15:23:06,967 | DEBUG | AstroPhotography |   key=Thumbnail JPEGInterchangeFormatLength  val=1832
2021-03-14 15:23:06,967 | DEBUG | AstroPhotography |   key=IFD2 ImageWidth                 val=384
2021-03-14 15:23:06,967 | DEBUG | AstroPhotography |   key=IFD2 ImageLength                val=256
2021-03-14 15:23:06,967 | DEBUG | AstroPhotography |   key=IFD2 BitsPerSample              val=[8, 8, 8]
2021-03-14 15:23:06,967 | DEBUG | AstroPhotography |   key=IFD2 Compression                val=JPEG (old-style)
2021-03-14 15:23:06,967 | DEBUG | AstroPhotography |   key=IFD2 PhotometricInterpretation  val=2
2021-03-14 15:23:06,967 | DEBUG | AstroPhotography |   key=IFD2 StripOffsets               val=308337
2021-03-14 15:23:06,968 | DEBUG | AstroPhotography |   key=IFD2 SamplesPerPixel            val=3
2021-03-14 15:23:06,968 | DEBUG | AstroPhotography |   key=IFD2 RowsPerStrip               val=256
2021-03-14 15:23:06,968 | DEBUG | AstroPhotography |   key=IFD2 StripByteCounts            val=294912
2021-03-14 15:23:06,968 | DEBUG | AstroPhotography |   key=IFD2 PlanarConfiguration        val=1
2021-03-14 15:23:06,968 | DEBUG | AstroPhotography |   key=IFD3 Compression                val=JPEG (old-style)
2021-03-14 15:23:06,968 | DEBUG | AstroPhotography |   key=IFD3 StripOffsets               val=603249
2021-03-14 15:23:06,968 | DEBUG | AstroPhotography |   key=IFD3 StripByteCounts            val=7721050
2021-03-14 15:23:06,968 | DEBUG | AstroPhotography |   key=EXIF ExposureTime               val=30
2021-03-14 15:23:06,969 | DEBUG | AstroPhotography |   key=EXIF FNumber                    val=28/5
2021-03-14 15:23:06,969 | DEBUG | AstroPhotography |   key=EXIF ExposureProgram            val=Manual
2021-03-14 15:23:06,969 | DEBUG | AstroPhotography |   key=EXIF ISOSpeedRatings            val=400
2021-03-14 15:23:06,969 | DEBUG | AstroPhotography |   key=EXIF ExifVersion                val=0221
2021-03-14 15:23:06,969 | DEBUG | AstroPhotography |   key=EXIF DateTimeOriginal           val=2014:01:18 21:48:57
2021-03-14 15:23:06,969 | DEBUG | AstroPhotography |   key=EXIF DateTimeDigitized          val=2014:01:18 21:48:57
2021-03-14 15:23:06,969 | DEBUG | AstroPhotography |   key=EXIF ComponentsConfiguration    val=
2021-03-14 15:23:06,970 | DEBUG | AstroPhotography |   key=EXIF ShutterSpeedValue          val=-321577/65536
2021-03-14 15:23:06,970 | DEBUG | AstroPhotography |   key=EXIF ApertureValue              val=162885/32768
2021-03-14 15:23:06,970 | DEBUG | AstroPhotography |   key=EXIF ExposureBiasValue          val=0
2021-03-14 15:23:06,970 | DEBUG | AstroPhotography |   key=EXIF MeteringMode               val=Average
2021-03-14 15:23:06,970 | DEBUG | AstroPhotography |   key=EXIF Flash                      val=Flash did not fire, compulsory flash mode
2021-03-14 15:23:06,970 | DEBUG | AstroPhotography |   key=EXIF FocalLength                val=30
2021-03-14 15:23:06,971 | DEBUG | AstroPhotography |   key=EXIF UserComment                val=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
2021-03-14 15:23:06,971 | DEBUG | AstroPhotography |   key=EXIF FlashPixVersion            val=0100
2021-03-14 15:23:06,971 | DEBUG | AstroPhotography |   key=EXIF ColorSpace                 val=Uncalibrated
2021-03-14 15:23:06,971 | DEBUG | AstroPhotography |   key=EXIF ExifImageWidth             val=3888
2021-03-14 15:23:06,971 | DEBUG | AstroPhotography |   key=EXIF ExifImageLength            val=2592
2021-03-14 15:23:06,971 | DEBUG | AstroPhotography |   key=Interoperability InteroperabilityIndex  val=R98
2021-03-14 15:23:06,971 | DEBUG | AstroPhotography |   key=Interoperability InteroperabilityVersion  val=[48, 49, 48, 48]
2021-03-14 15:23:06,972 | DEBUG | AstroPhotography |   key=EXIF InteroperabilityOffset     val=77610
2021-03-14 15:23:06,972 | DEBUG | AstroPhotography |   key=EXIF FocalPlaneXResolution      val=3888000/877
2021-03-14 15:23:06,972 | DEBUG | AstroPhotography |   key=EXIF FocalPlaneYResolution      val=432000/97
2021-03-14 15:23:06,972 | DEBUG | AstroPhotography |   key=EXIF FocalPlaneResolutionUnit   val=2
2021-03-14 15:23:06,972 | DEBUG | AstroPhotography |   key=EXIF CustomRendered             val=Normal
2021-03-14 15:23:06,972 | DEBUG | AstroPhotography |   key=EXIF ExposureMode               val=Manual Exposure
2021-03-14 15:23:06,972 | DEBUG | AstroPhotography |   key=EXIF WhiteBalance               val=Auto
2021-03-14 15:23:06,973 | DEBUG | AstroPhotography |   key=EXIF SceneCaptureType           val=Standard
2021-03-14 15:23:06,973 | DEBUG | AstroPhotography |   key=MakerNote FlashInfo             val=[100, 0, 0, 0]
2021-03-14 15:23:06,973 | DEBUG | AstroPhotography |   key=MakerNote ImageType             val=Canon EOS DIGITAL REBEL XTi
2021-03-14 15:23:06,973 | DEBUG | AstroPhotography |   key=MakerNote FirmwareVersion       val=Firmware 1.1.0
2021-03-14 15:23:06,973 | DEBUG | AstroPhotography |   key=MakerNote OwnerName             val=unknown
2021-03-14 15:23:06,973 | DEBUG | AstroPhotography |   key=MakerNote SerialNumber          val=1521037688
2021-03-14 15:23:06,973 | DEBUG | AstroPhotography |   key=MakerNote ModelID               val=EOS Digital Rebel XTi / 400D / Kiss Digital X
2021-03-14 15:23:06,973 | DEBUG | AstroPhotography |   key=MakerNote ThumbnailImageValidArea  val=[0, 159, 7, 112]
2021-03-14 15:23:06,974 | DEBUG | AstroPhotography |   key=MakerNote SerialNumberFormat    val=Format 2
2021-03-14 15:23:06,974 | DEBUG | AstroPhotography |   key=MakerNote LensModel             val=EF-S18-55mm f/3.5-5.6
2021-03-14 15:23:06,974 | DEBUG | AstroPhotography |   key=MakerNote InternalSerialNumber   val=H2485064
2021-03-14 15:23:06,974 | DEBUG | AstroPhotography |   key=MakerNote DustRemovalData       val=[]
2021-03-14 15:23:06,974 | DEBUG | AstroPhotography |   key=MakerNote CropInfo              val=[0, 0, 0, 0]
2021-03-14 15:23:06,974 | DEBUG | AstroPhotography |   key=MakerNote ColorSpace            val=Adobe RGB
2021-03-14 15:23:06,974 | DEBUG | AstroPhotography |   key=MakerNote Macromode             val=Normal
2021-03-14 15:23:06,975 | DEBUG | AstroPhotography |   key=MakerNote SelfTimer             val=0
2021-03-14 15:23:06,975 | DEBUG | AstroPhotography |   key=MakerNote Quality               val=Unknown
2021-03-14 15:23:06,975 | DEBUG | AstroPhotography |   key=MakerNote FlashMode             val=Flash Not Fired
2021-03-14 15:23:06,975 | DEBUG | AstroPhotography |   key=MakerNote ContinuousDriveMode   val=Single Or Timer
2021-03-14 15:23:06,975 | DEBUG | AstroPhotography |   key=MakerNote Unknown               val=0
2021-03-14 15:23:06,975 | DEBUG | AstroPhotography |   key=MakerNote FocusMode             val=Manual
2021-03-14 15:23:06,975 | DEBUG | AstroPhotography |   key=MakerNote RecordMode            val=CR2
2021-03-14 15:23:06,976 | DEBUG | AstroPhotography |   key=MakerNote ImageSize             val=Unknown
2021-03-14 15:23:06,976 | DEBUG | AstroPhotography |   key=MakerNote EasyShootingMode      val=Manual
2021-03-14 15:23:06,976 | DEBUG | AstroPhotography |   key=MakerNote DigitalZoom           val=None
2021-03-14 15:23:06,976 | DEBUG | AstroPhotography |   key=MakerNote Contrast              val=Normal
2021-03-14 15:23:06,976 | DEBUG | AstroPhotography |   key=MakerNote Saturation            val=Normal
2021-03-14 15:23:06,976 | DEBUG | AstroPhotography |   key=MakerNote Sharpness             val=Unknown
2021-03-14 15:23:06,976 | DEBUG | AstroPhotography |   key=MakerNote ISO                   val=Unknown
2021-03-14 15:23:06,976 | DEBUG | AstroPhotography |   key=MakerNote MeteringMode          val=Center-weighted
2021-03-14 15:23:06,977 | DEBUG | AstroPhotography |   key=MakerNote FocusType             val=Unknown
2021-03-14 15:23:06,977 | DEBUG | AstroPhotography |   key=MakerNote AFPointSelected       val=Unknown
2021-03-14 15:23:06,977 | DEBUG | AstroPhotography |   key=MakerNote ExposureMode          val=Manual
2021-03-14 15:23:06,977 | DEBUG | AstroPhotography |   key=MakerNote LensType              val=45
2021-03-14 15:23:06,977 | DEBUG | AstroPhotography |   key=MakerNote LongFocalLengthOfLensInFocalUnits  val=55
2021-03-14 15:23:06,977 | DEBUG | AstroPhotography |   key=MakerNote ShortFocalLengthOfLensInFocalUnits  val=18
2021-03-14 15:23:06,977 | DEBUG | AstroPhotography |   key=MakerNote FocalUnitsPerMM       val=1
2021-03-14 15:23:06,978 | DEBUG | AstroPhotography |   key=MakerNote FlashActivity         val=Did Not Fire
2021-03-14 15:23:06,978 | DEBUG | AstroPhotography |   key=MakerNote FlashDetails          val=Manual
2021-03-14 15:23:06,978 | DEBUG | AstroPhotography |   key=MakerNote AESetting             val=Unknown
2021-03-14 15:23:06,978 | DEBUG | AstroPhotography |   key=MakerNote ImageStabilization    val=Unknown
2021-03-14 15:23:06,978 | DEBUG | AstroPhotography |   key=MakerNote SpotMeteringMode      val=Unknown
2021-03-14 15:23:06,978 | DEBUG | AstroPhotography |   key=MakerNote ManualFlashOutput     val=n/a
2021-03-14 15:23:06,978 | DEBUG | AstroPhotography |   key=MakerNote FocalType             val=Unknown
2021-03-14 15:23:06,979 | DEBUG | AstroPhotography |   key=MakerNote FocalLength           val=907
2021-03-14 15:23:06,979 | DEBUG | AstroPhotography |   key=MakerNote WhiteBalance          val=Auto
2021-03-14 15:23:06,979 | DEBUG | AstroPhotography |   key=MakerNote SlowShutter           val=None
2021-03-14 15:23:06,979 | DEBUG | AstroPhotography |   key=MakerNote SequenceNumber        val=0
2021-03-14 15:23:06,979 | DEBUG | AstroPhotography |   key=MakerNote AFPointUsed           val=0
2021-03-14 15:23:06,979 | DEBUG | AstroPhotography |   key=MakerNote FlashBias             val=0 EV
2021-03-14 15:23:06,979 | DEBUG | AstroPhotography |   key=MakerNote SubjectDistance       val=222
2021-03-14 15:23:06,979 | DEBUG | AstroPhotography |   key=MakerNote FileNumber            val=36906
2021-03-14 15:23:06,980 | DEBUG | AstroPhotography |   key=MakerNote BracketMode           val=Off
2021-03-14 15:23:06,980 | DEBUG | AstroPhotography |   key=MakerNote BracketValue          val=0
2021-03-14 15:23:06,980 | DEBUG | AstroPhotography |   key=MakerNote BracketShotNumber     val=0
2021-03-14 15:23:06,980 | DEBUG | AstroPhotography |   key=MakerNote RawJpgQuality         val=n/a
2021-03-14 15:23:06,980 | DEBUG | AstroPhotography |   key=MakerNote RawJpgSize            val=Unknown
2021-03-14 15:23:06,981 | DEBUG | AstroPhotography |   key=MakerNote LongExposureNoiseReduction2  val=Off
2021-03-14 15:23:06,981 | DEBUG | AstroPhotography |   key=MakerNote WBBracketMode         val=Off
2021-03-14 15:23:06,981 | DEBUG | AstroPhotography |   key=MakerNote WBBracketValueAB      val=0
2021-03-14 15:23:06,981 | DEBUG | AstroPhotography |   key=MakerNote WBBracketValueGM      val=0
2021-03-14 15:23:06,981 | DEBUG | AstroPhotography |   key=MakerNote FilterEffect          val=None
2021-03-14 15:23:06,982 | DEBUG | AstroPhotography |   key=MakerNote ToningEffect          val=None
2021-03-14 15:23:06,982 | DEBUG | AstroPhotography |   key=MakerNote MacroMagnification    val=119
```
