# RdNapTrans

Python 3 implementation of RDNAPTRANS&trade;  
Converts Spatial Coordinates from RD_New (EPSG:28992 https://spatialreference.org/ref/epsg/amersfoort-rd-new/) to ETRS89 (EPSG:4258 https://spatialreference.org/ref/epsg/4258/) and back.  
Ported from the C version of RdNapTrans:  
See: https://zakelijk.kadaster.nl/transformatie-van-coordinaten  
The three correction grid files (nlgeo04.grd, x2c.grd, y2c.grd) needed for the grid interpolation (see https://nl.wikipedia.org/wiki/Rijksdriehoekscoördinaten) are included as base64 encoded strings. 
