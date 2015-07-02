rm *.mp4
ffmpeg -i log10_density_%d.png wd_nova_log10_density.mp4
ffmpeg -i log10_temperature_%d.png wd_nova_log10_temperature.mp4
ffmpeg -i x_velocity_%d.png wd_nova_x_velocity.mp4
ffmpeg -i y_velocity_%d.png wd_nova_y_velocity.mp4
