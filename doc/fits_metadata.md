## FITS Metadata

### FITS processing with the AstroPhotography package

The AstroPhotography package follows the standard astronomical data 
processing philosophy: 

- a series of separate tools representing sequential stages 
- each reading or writing one or more file-like entities 
- the file-like entities consist of array-like data with an associated
  set of metadata
- the metadata informs each stage of processing and is amended or added
  to by each stage
- in most cases missing metadata will cause processing to fail, or
  work less well. 
  
The purpose of the `ApAddMetadata` class (`ap_add_metadata` command line 
script) is to add to or update the metadata to the standard expected by
later `AstroPhotography` processing stages.

#### ApAddMetadata / ap_add_metadata

##### Metadata provided in iTelescope raw FITS files

The FITS headers of FITS images obtained from iTelescope user observsing
sessions typically contain the observing time, some fixed telescope 
characteristics, and the detector temperature

For example, the file `raw-T31-davestrickland-webb-20211229-024225-Green-BIN1-W-300-001.fit` 
contains the following:

```
# HDU 0 in raw-T31-davestrickland-webb-20211229-024225-Green-BIN1-W-300-001.fit:
SIMPLE  =                    T                                                  
BITPIX  =                   16 /8 unsigned int, 16 & 32 int, -32 & -64 real     
NAXIS   =                    2 /number of axes                                  
NAXIS1  =                 3056 /fastest changing axis                           
NAXIS2  =                 3056 /next to fastest changing axis                   
BSCALE  =   1.0000000000000000 /physical = BZERO + BSCALE*array_value           
BZERO   =   32768.000000000000 /physical = BZERO + BSCALE*array_value           
DATE-OBS= '2021-12-29T15:42:46.53' /YYYY-MM-DDThh:mm:ss observation, UT         
EXPTIME =   300.00000000000000 /Exposure time in seconds                        
EXPOSURE=   300.00000000000000 /Exposure time in seconds                        
SET-TEMP=  -25.000000000000000 /CCD temperature setpoint in C                   
CCD-TEMP=  -25.000000000000000 /CCD temperature at start of exposure in C       
XPIXSZ  =   12.000000000000000 /Pixel Width in microns (after binning)          
YPIXSZ  =   12.000000000000000 /Pixel Height in microns (after binning)         
XBINNING=                    1 /Binning factor in width                         
YBINNING=                    1 /Binning factor in height                        
XORGSUBF=                    0 /Subframe X position in binned pixels            
YORGSUBF=                    0 /Subframe Y position in binned pixels            
READOUTM= '8 MHz (RBI Flood)' / Readout mode of image                           
FILTER  = 'Green   ' /          Filter used when taking image                   
IMAGETYP= 'Light Frame' /       Type of image                                   
FOCALLEN=   2259.0000000000000 /Focal length of telescope in mm                 
APTDIA  =   510.00000000000000 /Aperture diameter of telescope in mm            
APTAREA =   204282.06798434258 /Aperture area of telescope in mm^2              
SBSTDVER= 'SBFITSEXT Version 1.0' /Version of SBFITSEXT standard in effect      
SWCREATE= 'MaxIm DL Version 6.27 211016 1NURE' /Name of software                
SWSERIAL= '1NURE-W5KP0-V31P2-TUEQ3-MM333-3P' /Software serial number            
SITELAT = '-31 16 24' /         Latitude of the imaging location                
SITELONG= '149 04 11' /         Longitude of the imaging location               
JD      =   2459578.1547052083 /Julian Date at time of exposure                 
JD-HELIO=   2459578.1556659853 /Heliocentric Julian Date at time of exposure    
OBJECT  = '        '                                                            
TELESCOP= '        ' /          telescope used to acquire this image            
INSTRUME= 'FLI     ' /          instrument or camera used                       
OBSERVER= '        '                                                            
NOTES   = '        '                                                            
ROWORDER= 'TOP-DOWN' /          Image write order, BOTTOM-UP or TOP-DOWN        
FLIPSTAT= '        '                                                            
CSTRETCH= 'Medium  ' /          Initial display stretch mode                    
CBLACK  =                 2538 /Initial display black level in ADUs             
CWHITE  =                 4291 /Initial display white level in ADUs             
PEDESTAL=                    0 /Correction to add for zero-based ADU            
SWOWNER = 'Brad Moore' /        Licensed owner of software       
```

Note that `OBJECT`, `TELESCOP`, and `OBSERVER` are blank in the FITS 
header, but are present in the file name. 

##### Metadata provided in iTelescope "Premium" FITS files

In comparison, the FITS headers of iTelescope premium images are 
impoverished and the file names no longer follow the normal user naming
conversion. 

For example, the file `Calibrated-iTelescope-Whale_Galaxy_NGC_55-900s-Red-5.fits` 
has the following FITS header:

```
# HDU 0 in Calibrated-iTelescope-Whale_Galaxy_NGC_55-900s-Red-5.fits:
SIMPLE  =                    T / conforms to FITS standard                      
BITPIX  =                  -32 / array data type                                
NAXIS   =                    2 / number of array dimensions                     
NAXIS1  =                 4096                                                  
NAXIS2  =                 4096                                                  
EXTEND  =                    T                                                  
OWNER   = 'iTelescope'                                                          
GAIN    =                 1.31 / gain or ISO depending on instrument            
EXPTIME =                900.0 / Exposure time in seconds                       
FILTER  = 'Red     '           / Filter type used                               
HIERARCH OBSTARGET = 'Whale Galaxy NGC 55' / Object imaged                      
SRCLINK = 'obqgw8i8'           / Original source file ID  
```
