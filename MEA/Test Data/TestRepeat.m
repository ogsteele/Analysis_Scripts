ElPos_SpikeNum = [(ElPos)',(SpikeNums)'];
    
Corrected_ElPos_SpikeNum = zeros(12,2);
Corrected_ElPos_SpikeNum(:,1) = [21,31,12,22,32,42,13,23,33,43,24,34]';
Empty = zeros(0,1);

ElInd = zeros(12,1);

    El21 = find(ElPos_SpikeNum(:,1) == 21);
        if El21 == Empty
            Corrected_ElPos_SpikeNum(1,2) = 0;
        else
            Corrected_ElPos_SpikeNum(1,2) = ElPos_SpikeNum(El21,2);
        end
    El31 = find(ElPos_SpikeNum(:,1) == 31);
        if El31 == Empty
            Corrected_ElPos_SpikeNum(2,2) = 0;
        else
            Corrected_ElPos_SpikeNum(2,2) = ElPos_SpikeNum(El31,2);
        end
    El12 = find(ElPos_SpikeNum(:,1) == 12);
        if El12 == Empty
            Corrected_ElPos_SpikeNum(3,2) = 0;
        else
            Corrected_ElPos_SpikeNum(3,2) = ElPos_SpikeNum(El12,2);
        end
    El22 = find(ElPos_SpikeNum(:,1) == 22);
        if El22 == Empty
            Corrected_ElPos_SpikeNum(4,2) = 0;
        else
            Corrected_ElPos_SpikeNum(4,2) = ElPos_SpikeNum(El22,2);
        end
    El32 = find(ElPos_SpikeNum(:,1) == 32);
        if El32 == Empty
            Corrected_ElPos_SpikeNum(5,2) = 0;
        else
            Corrected_ElPos_SpikeNum(5,2) = ElPos_SpikeNum(El32,2);
        end
    El42 = find(ElPos_SpikeNum(:,1) == 42);
        if El42 == Empty
            Corrected_ElPos_SpikeNum(6,2) = 0;
        else
            Corrected_ElPos_SpikeNum(6,2) = ElPos_SpikeNum(El42,2);
        end
    El13 = find(ElPos_SpikeNum(:,1) == 13);
        if El13 == Empty
            Corrected_ElPos_SpikeNum(7,2) = 0;
        else
            Corrected_ElPos_SpikeNum(7,2) = ElPos_SpikeNum(El13,2);
        end
    El23 = find(ElPos_SpikeNum(:,1) == 23);
        if El23 == Empty
            Corrected_ElPos_SpikeNum(8,2) = 0;
        else
            Corrected_ElPos_SpikeNum(8,2) = ElPos_SpikeNum(El23,2);
        end
    El33 = find(ElPos_SpikeNum(:,1) == 33);
        if El33 == Empty
            Corrected_ElPos_SpikeNum(9,2) = 0;
        else
            Corrected_ElPos_SpikeNum(9,2) = ElPos_SpikeNum(El33,2);
        end
    El43 = find(ElPos_SpikeNum(:,1) == 43);
        if El43 == Empty
            Corrected_ElPos_SpikeNum(10,2) = 0;
        else
            Corrected_ElPos_SpikeNum(10,2) = ElPos_SpikeNum(El43,2);
        end
    El24 = find(ElPos_SpikeNum(:,1) == 24);
        if El24 == Empty
            Corrected_ElPos_SpikeNum(11,2) = 0;
        else
            Corrected_ElPos_SpikeNum(11,2) = ElPos_SpikeNum(El24,2);
        end
    El34 = find(ElPos_SpikeNum(:,1) == 34);
        if El34 == Empty
            Corrected_ElPos_SpikeNum(12,2) = 0;
        else
            Corrected_ElPos_SpikeNum(12,2) = ElPos_SpikeNum(El34,2);
        end
        
        
        
%for x = 1:12
    %ElInd(x) = find(ElPos_SpikeNum(:,1) == (Corrected_ElPos_SpikeNum(x,1)));
       % if ElInd(x) == Empty
           % ElInd(x) = 0;
           % Corrected_ElPos_SpikeNum(x,2) = 0;
       % else
           % Corrected_ElPos_SpikeNum(x,2) = ElPos_SpikeNum(ElInd(x),2);
       % end
% end