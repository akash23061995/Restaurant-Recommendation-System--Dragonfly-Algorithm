clc;
clear;
close all;

possibility=7;
while(1)
    xxx=menu('RMS','Load Data','Preprocess','Apply Rough Set','Profile Aggregation','Apply DA','Exit');
    switch(xxx)
        case 1
            loaddata
        case 2
            preprocess
        case 3
            aprioritest
        case 4
            PAA
        case 5
            applyda
        case 6
            break;
    end
end

