rm *.mp4
ffmpeg -i log10_density_%d.png -codec:v libvpx wd_nova_log10_density.webm
ffmpeg -i log10_temperature_%d.png -codec:v libvpx wd_nova_log10_temperature.webm
ffmpeg -i x_velocity_%d.png -codec:v libvpx wd_nova_x_velocity.webm
ffmpeg -i y_velocity_%d.png -codec:v libvpx wd_nova_y_velocity.webm
