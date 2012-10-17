REM Testing Unwrap-wraproutines

release\Panorama sphrWarp images\pano1_0008.tga images\t_warp1.tga 595 -0.15 0.0

REM Test align
release\Panorama alignPair images\warp08.key images\warp09.key images\match-08-09.txt 200 4 sift
release\Panorama alignPairHomography images\warp08.key images\warp09.key images\match-08-09.txt 200 4 sift