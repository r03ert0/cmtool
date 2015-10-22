cmtool
======

<p>cmtool is a command line tool written in standard C to query the
Brain Coactivation Map data
(<a href="http://www.nitrc.org/projects/cmap/">available from NITRC</a>).</p>
<p>Alternatively, you can also use our interactive viewer <a href="https://github.com/r03ert0/CoactivationMap.app/">CoactivationMap</a>.</p>

## Arguments
````
    argv[1]: cmapdir                          | cmapdir=coactivation map directory (contains
                                              | defaults.txt, sum.img, brainspell.xml)
    argv[2,3]:{mni x,y,z|vol i,j,k|roi path|  | mni x,y,z= seed coordinates in MNI space,
                tag name|vol6 i,j,k,x,y,z|    | vol i,j,k= seed coordinates in volume index
                                              | space,
                mnifile path}                 | roi {average} path= path to ROI file to use
                                              | as multiple seeds or average through all seeds
                                              | tag name= name of the tag to pull
                                              | vol6 i,j,k are seed coordinates, x,y,z are
                                              | space coordinates, 3 coordinates must have a
                                              | number, the others will be scanned
                                              | mnifile path= path to text file containing
                                              | seeds in MNI space,
    argv[4...] commands:
            -swap                             | swap endianness in coactivation map files
                                              | (coincidence, sum)
            -average                          | average over the ROI
            -convert {lr|mi|logp|t|phir|dct}  | convert map to lr=likelihood ratio, mi=mutual
                                              | information, logp=log p-value of likelihood
                                              | ratio test, t=t-value, phir=phi-correlation,
                                              |     dct=discrete cosinus transform
            -out root                         | root=root name for output files (end with /
                                              | for directory)
            -meshdisplay {name|code}          | display mesh tags as names or codes (default
                                              | is "codes")
            -meshfilter    root               | in article selection, only report mesh heading
                                              | codes (not names) under root (default "F")
            -coverseed                        | in article selection, only report articles that
                                              | cover the seed node (default, NO)
            -coverage thrs                    | in article selection, minimum coverage (from 0
                                              | to 100) of network nodes to report an article
                                              | (default, 0)
            -save map                         | save the map
        [DEL]    -save average                | save the average map (only for argv[2]=roi)
            -save peaks thrs R                | save list of map peaks detected with threshold
                                              | thrs and radius R
            -save articles thrs R             | save list of articles intersecting the map
                                              | peaks detected with threshold thrs and radius R
            -save topmesh                     | rank mesh tags based on a map-weighted sum
                                              | throughout all papers. Save the top 100
            -saveroiarticles                  | save list of articles intersecting the ROI
                                              | (arg[2,3]=roi path)
````

## Examples

Get the map corresponding to voxel indices 6,24,23 as a likelihood ratio volume
````
./cmtool coincidences2013 -vol 6,24,23 -convert lr -out ~/Desktop/ -save map
````

Get the map corresponding to voxel 6,24,23 as a phi-correlation volume
````
./cmtool coincidences2013 -vol 6,24,23 -convert phir -out ~/Desktop/ -save map
````

Compute all the maps for voxels in the left insula region of interest, convert them to phi-correlations, average them and save
````
./cmtool coincidences2013 -roi 01roi/left_insula_ROI.hdr -convert phir -average -out ~/Desktop/ -save map
````

Compute the phi-correlation map corresponding to voxel 6,24,23 and save the peaks over 0.2 within a neighbourhood of size 4 voxels
````
./cmtool coincidences2013 -vol 6,24,23 -convert phir -out ~/Desktop/ -save peaks 0.2 4
````

Compute the phi-correlation map corresponding to voxel 6,24,23, compute the peaks over 0.2 within a neighbourhood of size 4 voxels, save the tags of the resulting articles
````
./cmtool coincidences2013 -vol 6,24,23 -convert phir -out ~/Desktop/ -save articles 0.2 4
````

Compute all the maps for voxels in the left insula region of interest, convert them to phi-correlations, average them and save the list of peaks over 0.1 within a neighbourhood of 4
````
./cmtool coincidences2013 -roi 01roi/left_insula_ROI.hdr -convert phir -average -out ~/Desktop/ -save peaks 0.1 4
````

More examples:

````
./cmtool coincidences2013 -roi 01roi/left_insula_ROI.hdr -convert phir -average -coverseed -out ~/Desktop/ -save articles 0.1 4
./cmtool coincidences2013 -roi 01roi/left_insula_ROI.hdr -convert phir -average -coverseed -coverage 50 -out ~/Desktop/ -save articles 0.1 4
./cmtool coincidences2013 -roi 01roi/left_insula_ROI.hdr -out ~/Desktop/ -saveroiarticles
./cmtool coincidences2013 -vol 6,24,23 -convert dct -out ~/Desktop/ -save map
./cmtool coincidences2013 -roimask.hdr -convert dct -out ~/Desktop/cmap-dct/ -average -save map
./cmtool coincidences2013 -vol 6,24,23 -convert phir -out ~/Desktop/ -save topmesh
./cmtool coincidences2013 -vol 6,24,23 -convert phir -out ~/Desktop/test -save articles 0.1 4
./cmtool coincidences2013 -vol 6,24,23 -convert phir -out ~/Desktop/test -meshfilter F03 -save articles 0.1 4
./cmtool coincidences2013 -vol 6,24,23 -convert phir -out ~/Desktop/test -meshfilter F03 -meshdisplay name -save articles 0.1 4
````
