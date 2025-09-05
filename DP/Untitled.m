theta1 = zeros(1,8001);
for i=1:8001
    if i < 122
        theta1(i) = 0;
    elseif i < 322
        theta1(i) = 0.07448;
    elseif i < 442
        theta1(i) = 0.7505;
    elseif i < 602
        theta1(i) = 0.8479;
    elseif i< 802
        theta1(i) = 1.066;
    elseif i < 882
        theta1(i) = 0.3667;
    elseif i < 1242
        theta1(i) = 0;
    elseif i < 1442
        theta1(i) = -0.3495;
    elseif i < 1682
        theta1(i) = 0;
    elseif i < 1882
        theta1(i) = -0.1375;
    elseif i < 2522
        theta1(i) = 0;
    elseif i < 2722
        theta1(i) = 0.3209;
    elseif i < 2882
        theta1(i) = -0.2177;
    elseif i < 3482
        theta1(i) = 0;
    elseif i < 3682
        theta1(i) = -0.2464;
    elseif i < 3842
        theta1(i) = 0.3094;
    elseif i < 4522
        theta1(i) = 0;
    elseif i < 4768
        theta1(i) = 0.3266;
    elseif i < 4842
        theta1(i) = 0;
    elseif i < 5083
        theta1(i) = 0.3352;
    elseif i < 5122
        theta1(i) = 0;
    elseif i < 5322
        theta1(i) = 0.5042;
    elseif i < 5522
        theta1(i) = 0;
    elseif i < 5722
        theta1(i) = 0.8078;
    elseif i < 5842
        theta1(i) = 0.03438;
    elseif i < 6722
        theta1(i) = 0;
    elseif i < 6922
        theta1(i) = -0.4813;
    elseif i < 7002
        theta1(i) = -0.2979;
    elseif i < 7442
        theta1(i) = 0;
    elseif i < 7642
        theta1(i) = 0.3953;
    else
        theta1(i) = 0;
    end
end