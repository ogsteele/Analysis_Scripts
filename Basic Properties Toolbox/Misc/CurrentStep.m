function [IR] = CurrentStep(filenameabf)
CS = abfload(filenameabf);
IR = (mean(CS(49112:51612,:,1) - CS(49112:51612,:,2)))*-100;
end

