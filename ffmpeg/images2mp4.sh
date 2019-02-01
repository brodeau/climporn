#!/bin/bash

# Laurent Brodeau, 2017-2019

# Output framerate for the video to generate in "fps":
#FRAMERATE_OUT="30000/1001" ; #NTSC (~29.97 fps)
FRAMERATE_OUT="25" ; # Normal ?

# Defaults:
FIMG="png"
HEIGHT="1080"
CODEC="x264"
FRAMERATE_IN=25
CRF=20
PRESET="medium"
FVID="mp4"
NTHRD="1"
DROPFRAMES=0
DIR_OUT="."

usage()
{
    echo
    echo "USAGE: `basename ${0}` -i <common_prefix_image_files> (options)"
    echo
    echo "   Available options are:"
    echo "      -t: format of images (default='${FIMG}')"
    echo "      -h: height (pixels) of video to create (default='${HEIGHT}')"
    echo "      -c: codec for video (default='${CODEC})' [x264, x265, ...]"
    echo "      -f: input framerate == images per second (default='${FRAMERATE_IN}')"
    echo "      -C: CRF value, 0 is lossless and 51 worse possible (default='${CRF}') [ffmpeg default=23]"
    echo "      -p: preset for encoding (default='${PRESET}') [fast, medium, slow, veryslow]"
    echo "      -v: video format (default='${FVID}') [mp4,webm,...]"
    echo "      -d: frequency for dropping frames (ex: -d 2 would speed up video by 2 by dropping every other frame)"
    echo "      -D: directory (path) in which to create the video if elsewhere than current directory"
    echo "      -n: number of threads (default='${NTHRD}')"
    echo
    exit
}


if [ "`ffmpeg -codecs 2>/dev/null | grep libx264`" = "" ]; then
    echo "Dude! Use a 'ffmpeg' that has 'libx264' support! Sorry..."
    exit
fi


while getopts i:t:h:c:f:C:p:v:d:D:n:h option; do
    case $option in
        i) FPREF=${OPTARG};;
        t) FIMG=${OPTARG};;
        h) HEIGHT=${OPTARG};;
        c) CODEC=${OPTARG};;
        f) FRAMERATE_IN=${OPTARG};;
        C) CRF=${OPTARG};;
        p) PRESET=${OPTARG};;
        v) FVID=${OPTARG};;
        d) DROPFRAMES=${OPTARG};;
        D) DIR_OUT=${OPTARG};;
        n) NTHRD=${OPTARG};;
        h)  usage ;;
        \?) usage ;;
    esac
done

if [ "${FPREF}" = "" ]; then usage; fi


# Video filter stuff:
VFLTR="-vf scale='-2:${HEIGHT}'"

if [ ${DROPFRAMES} -gt 1 ]; then
    if [ "${FRAMERATE_IN}" != "${FRAMERATE_OUT}" ]; then
        echo "PROBLEM: chose either option!"
        echo "   => if you drop frames with '-d n' then ensure that FRAMERATE_IN == FRAMERATE_OUT !"
        exit
    fi    
    if [ ${DROPFRAMES} -eq 2 ]; then
        VFLTR="${VFLTR},setpts='0.5*PTS'" ; # drop every other frame => speed of video x 2
    elif [ ${DROPFRAMES} -eq 3 ]; then
        VFLTR="${VFLTR},setpts='0.3*PTS'"
    elif [ ${DROPFRAMES} -eq 4 ]; then
        VFLTR="${VFLTR},setpts='0.25*PTS'"
    elif [ ${DROPFRAMES} -eq 5 ]; then
        VFLTR="${VFLTR},setpts='0.2*PTS'"
    elif [ ${DROPFRAMES} -eq 10 ]; then
        VFLTR="${VFLTR},setpts='0.1*PTS'"
    else
        echo "ERROR: we dont do your DROPFRAMES = ${DROPFRAMES} !"; exit
    fi
fi


#Codec stuff:
if [ "${FVID}" = "mp4" ]; then
    VC="-c:v libx264 -profile:v high444"
    info="${CODEC}_${HEIGHT}px"
elif [ "${FVID}" = "webm" ]; then
    VC="-c:v libvpx"
    #VC="-c:v libvpx -pix_fmt yuva420p -metadata:s:v:0 alpha_mode=\"1\""
    info="vpx_${HEIGHT}px"
elif [ "${FVID}" = "mov" ]; then
    VC="-c:v h264_nvenc"
    info="h264_${HEIGHT}px"
else
    echo "Boo!" ; exit
fi


fo="${DIR_OUT}/movie_${FPREF}_${info}_${FRAMERATE_IN}fps_crf${CRF}.${FVID}"

rm -f ${fo}

echo
echo " fo = ${fo} !!!"
#exit

echo
echo "ffmpeg -f image2 -threads ${NTHRD} \
-pattern_type glob -r ${FRAMERATE_IN} -i ${FPREF}*.${FIMG} \
${VC} -preset ${PRESET} \
-crf ${CRF} -refs 16 ${VFLTR} \
-pix_fmt yuv420p \
-r ${FRAMERATE_OUT} ${fo}"
echo


#ffmpeg -f image2 -threads ${NTHRD} -framerate ${FRAMERATE_IN} -r ${FRAMERATE_OUT} \
#       -pattern_type glob -i "${FPREF}*.${FIMG}" \
#       ${VC} -preset ${PRESET} \
#       -crf ${CRF} -refs 16 ${VFLTR} \
#       -pix_fmt yuv420p \
#       ${fo}

ffmpeg -f image2 -threads ${NTHRD} \
       -pattern_type glob -r ${FRAMERATE_IN} -i "${FPREF}*.${FIMG}" \
       ${VC} -preset ${PRESET} \
       -crf ${CRF} -refs 16 ${VFLTR} \
       -pix_fmt yuv420p \
       -r ${FRAMERATE_OUT} ${fo}

echo
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
