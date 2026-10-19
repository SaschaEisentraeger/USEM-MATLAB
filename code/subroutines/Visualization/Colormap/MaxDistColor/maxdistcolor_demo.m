%%% MAXDISTCOLOR Demo Script %%%
% Plot the output colormap in the UCS, as a distance plot, and as colorbars with colornames.
% Some examples use CAM02 colorspace functions, which must be downloaded separately.
% The CIELAB colorspace can be used, but the results are suboptimal.
%fun = @srgb_to_Lab;  % CIELAB = not recommended colorspace.
fun = @(m)Lab_to_DIN99(srgb_to_Lab(m)); % DIN99 = good colorspace (close colors).
%fun = @srgb_to_Lab; % CAM02-UCS = best colorspace (close colors).
%fun = @(m)srgb_to_Lab(m,true,'LCD'); % CAM02-LCD = best colorspace (distant colors).
%%% Generate colormap of distinctive colors:
[rgb,ucs] = maxdistcolor(9,fun);
%[rgb,ucs] = maxdistcolor(9,fun,'exc',[]);
%[rgb,ucs] = maxdistcolor(9,fun,'Cmin',0.5,'Cmax',0.6);
%[rgb,ucs] = maxdistcolor(9,fun,'Lmin',0.4,'Lmax',0.6);
%[rgb,ucs] = maxdistcolor(9,fun,'inc',[0,0,0;1,0,1],'exc',[0,1,0]);
%[rgb,ucs] = maxdistcolor(9,fun,'sort','longest','disp','verbose');
%[rgb,ucs] = maxdistcolor(9,fun, 'bitR',8,'bitG',8,'bitB',8); % 24bit RGB -> slow!
%[rgb,ucs] = maxdistcolor(63,fun,'bitR',2,'bitG',2,'bitB',2); % default excludes white.
%[rgb,ucs] = maxdistcolor(64,fun,'bitR',2,'bitG',2,'bitB',2, 'exc',[]); % entire RGB gamut.
N = size(rgb,1);
%%% Get colorspace axes:
f2s = [regexp(func2str(fun),'(?<=to_)([JL]ab|DIN99)','match','once'),'XXX'];
tmp = struct('Lab',{'J''','a''','b'''},'Lab',{'L*','a*','b*'},'DIN',{'L_{99}','a_{99}','b_{99}'},'XXX',{'L?','a?','b?'});
xyz = {tmp.(f2s(1:3))};
f2s = regexprep(f2s,{'^XXX$','XXX'},{'???',''});
%%% Plot color distance matrix:
figure();
for k = 1:N
	dst = sqrt(sum(bsxfun(@minus,ucs,ucs(k,:)).^2,2));
	scatter3(k*ones(1,N),1:N,dst, 123, rgb,...
		'MarkerFaceColor',rgb(k,:), 'LineWidth',2.8, 'Marker','o')
	hold on
end
title(sprintf('Euclidean Distances in %s Colorspace',f2s))
zlabel('Euclidean Distance')
ylabel('Colormap Index')
xlabel('Colormap Index')
set(gca,'XTick',1:N,'YTick',1:N)
%%% Plot colors in UCS:
figure();
scatter3(ucs(:,3),ucs(:,2),ucs(:,1), 256, rgb, 'filled')
text(ucs(:,3),ucs(:,2),ucs(:,1),cellstr(num2str((1:N).')), 'HorizontalAlignment','center')
%%% Plot outline of RGB cube:
M = 23;
[X,Y,Z] = ndgrid(linspace(0,1,M),0:1,0:1);
mat = fun([X(:),Y(:),Z(:);Y(:),Z(:),X(:);Z(:),X(:),Y(:)]);
X = reshape(mat(:,3),M,[]);
Y = reshape(mat(:,2),M,[]);
Z = reshape(mat(:,1),M,[]);
line(X,Y,Z,'Color','k')
axis('equal')
title(sprintf('Colormap in %s Colorspace',f2s))
zlabel(sprintf('Lightness (%s)',xyz{1}))
ylabel(sprintf('Hue Dim.1 (%s)',xyz{2}))
xlabel(sprintf('Hue Dim.2 (%s)',xyz{3}))
%%% Plot colorband image:
figure()
image(permute(rgb,[1,3,2]))
title('Colormap in Colorbands')
ylabel('Colormap Index')
set(gca,'XTick',[], 'YTick',1:N, 'YDir','normal')
%%% Add colornames (if COLORNAMES is available):
try %#ok<TRYNC>
	text(ones(1,N), 1:N, colornames('SVG',rgb),...
		'HorizontalAlignment','center', 'BackgroundColor','white')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%maxdistcolor_demo