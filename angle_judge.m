function flag = angle_judge(Th, Th_min, Th_max)

if Th >= Th_min && Th <= Th_max
    flag = 1;
else
    flag = 0;
end

end