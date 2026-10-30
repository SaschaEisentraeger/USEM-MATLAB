function ProgressBarUpdate(h,name,i,NumTimeStep,nUpdate)
%PROGRESSBARUPDATE Redraw the waitbar of ProgressBar after time step i.
%
%   ProgressBarUpdate(h, name, i, NumTimeStep, nUpdate) updates the waitbar
%   h every nUpdate steps (ProgressBar calls it through pb.update(i)).

if mod(i,nUpdate) == 0
    waitbar(i/NumTimeStep,h,[name ': ' ...
            num2str(i/NumTimeStep*100,'%5.1f') '%'])
end

end
