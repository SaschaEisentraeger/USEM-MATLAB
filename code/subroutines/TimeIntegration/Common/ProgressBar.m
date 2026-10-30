function pb = ProgressBar(name,NumTimeStep,enabled)
%PROGRESSBAR Waitbar of a time integration loop.
%
%   pb = ProgressBar(name, NumTimeStep) opens a waitbar titled name and
%   returns a struct with the function handles
%     pb.update(i)  show the progress after time step i (the bar is
%                   redrawn at most 1000 times to keep the loop fast)
%     pb.close()    close the waitbar
%   pb = ProgressBar(name, NumTimeStep, false) returns handles that do
%   nothing.

if nargin < 3
    enabled = true;
end
if ~enabled
    pb.update = @(i) [];
    pb.close  = @() [];
    return
end

h = waitbar(0,[name ': ']);
nUpdate = ceil(NumTimeStep/1000);
pb.update = @(i) ProgressBarUpdate(h,name,i,NumTimeStep,nUpdate);
pb.close  = @() close(h);

end
