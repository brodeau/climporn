#!/bin/bash

#PRESET="veryslow"
#PRESET="slow"
PRESET="medium"
#PRESET="ultrafast"

FPS=12

#CRF=20
CRF=18
#CRF=0

TYPE="mp4"
#TYPE="webm"
#TYPE="mov"


if [ "${5}" = "" ]; then
    echo "USAGE: `basename $0` <begining_files> <format (jpg,png,...)> <height video (pixels)> <fps> <crf>"
    exit
fi

SCALE="-vf scale='-2:${3}'"
FPS=${4}
CRF=${5}

#Codec stuff

if [ "${TYPE}" = "mp4" ]; then
    VC="-c:v libx264 -profile:v high444"
    info="x264_${3}px"
elif [ "${TYPE}" = "webm" ]; then
    VC="-c:v libvpx"
    #VC="-c:v libvpx -pix_fmt yuva420p -metadata:s:v:0 alpha_mode=\"1\""
    info="vpx_${3}px"
elif [ "${TYPE}" = "mov" ]; then
    VC="-c:v h264_nvenc"
    info="h264_${3}px"
else
    echo "Boo!" ; exit
fi



fo="movie_${1}_${info}_${FPS}fps_crf${CRF}.${TYPE}"

rm -f ${fo}



ffmpeg -f image2 -framerate ${FPS} \
       -pattern_type glob -i "${1}*.${2}" \
       ${VC} -preset ${PRESET} \
       -crf ${CRF} -refs 16 ${SCALE} \
       -pix_fmt yuv420p \
       ${fo}


echo
echo " *** Check ${fo} !!!"
echo
exit 0


#ffmpeg -f image2 -r 30 -i %09d.jpg -vcodec libx264 -profile:v high444 -refs 16 -crf 0 -preset ultrafast a.mp4

#Explanation of options:

#    -f image2 - tells ffmpeg to select a group of images
#    -r 30 - tells ffmpeg to encode at 30 frames (or images) per second (change this to your desired framerate)
#    -vcodec libx264 - tells ffmpeg to output to a H.264 compliant file
#    -profile:v high444 - tells libx264 to use the High 4:4:4 Predictive Lossless profile, allowing lossless encoding
#    -refs 16 - tells libx264 to have 16 images stored in a buffer, so that they may be referenced by other images in the video
#    -crf 0 - tells libx264 to perform a lossless encode
#    -preset ultrafast - tells libx264 to prioritise encoding speed over output file size
