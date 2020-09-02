ffmpeg -framerate 15 -i %02d.jpeg -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p lorenz.mp4

 
