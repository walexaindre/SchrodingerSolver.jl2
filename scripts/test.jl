using ProgressMeter
@showprogress "Comp" showspeed=true barlen=30 for j in 1:100
    @showprogress "V" offset=1 showspeed=true barlen=30 for i in 1:10
        sleep(0.25)
    end
end