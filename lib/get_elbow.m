% -------------------------------------------------------------------------
% GET_ELBOW calculates the location of the elbow for a given vector of
% values by taking the min or max distance from the origin depending on if
% the elbow is convex or concave, respectively. If the vector is neither,
% the middle index is given.
% -------------------------------------------------------------------------

function thresh = get_elbow(v)

% Check for all zeros weights.
if sum(v) == 0
    thresh = 0;
    return
end

% Scale vectors between 0 and 1.
y = (v - min(v))/(max(v) - min(v));
x = 0:(length(v) - 1);
x = x/max(x);

% Calculate distance from the origin.
d = sqrt(x.^2 + y.^2);

% Calculate number of points above and below unit circle.
above = sum(d > 1);
below = sum(d < 1);

% Determine threshold.
if above > below
    [~, thresh] = max(d);
elseif above < below
    [~, thresh] = min(d);
else
    thresh = length(d)/2;
end

% Check that threshold is not a zero weighted edge.
if v(thresh) == 0
    thresh = thresh - 1;
end

% Optional plotting.
% h = figure;
% h.Position = [300 300 1200 600];
% 
% subplot(1,2,1);
% scatter(x,y,20,'.k');
% axis square;
% hold on;
% plot(x,y,'k');
% plot(x,sqrt(1 - x.^2),'k:');
% scatter(x(thresh),y(thresh),25,'filled','r');
% plot([0 x(thresh)],[0 y(thresh)],'r');
% 
% subplot(1,2,2);
% i = 1:length(d);
% scatter(i,d,20,'.k');
% axis square;
% hold on;
% plot(i,d,'k');
% scatter(i(d > 1),d(d > 1),20,'.r');
% scatter(i(d < 1),d(d < 1),20,'.b');
% plot([0 length(d)],[1 1],':','Color',[0.5 0.5 0.5]);
% plot(thresh, d(thresh),'gd','MarkerSize',3,'MarkerFaceColor','g');
% plot([0 thresh thresh],[d(thresh) d(thresh) 0],'g:');
% axis([0 length(d) min(d) - 0.01 max(d) + 0.01])
% 
% close all

end