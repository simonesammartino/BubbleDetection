% Script to detect bubbling on water surface, and define number of bubbles
% and their size, from a video footage of the methane activity of a small
% pond. This version of the code is prepared to be run on a short test
% video footage provided along with the code.

% Author: Simone Sammartino
% Physical Oceanography Group
% University of Málaga
% email: ssammartino@ctima.uma.es
% -----------------------------------------------------------------------%

clear
close all
clc

lim = [1790 2540 850 1270];
nj = lim(2)-lim(1)+1;
ni = lim(4)-lim(3)+1;
nfm = 70; % Maximum number of frames analyzed
res = 2.1^2; % Ground resolution (GSD) mm/px
msx = 20; % Max bubble diameter
buf = round(msx/res); % Max number of pixels in neighboors

% Video file
vid = VideoReader("Test.mp4");

%% Graphical interface
figure('Position',[5 20 1750 1000],'Color','w')

% Image
axes('Position',[0.01 0.02 0.84 0.8])
fra = read(vid,[1 1]);
fra = fra(lim(3):lim(4),lim(1):lim(2),1);
img = fra;
pfr = imshow(img);
hold on
plt = scatter(nan,nan,'oy','Linewidth',2);

% Histogram
axes('Position',[0.882 0.07 0.11 0.75])
hi = 0.5:1:20.5;
brh = barh(hi(1:end-1)+0.5*diff(hi(1:2)),ones(size(hi,1),size(hi,2)-1));
grid on
axis([0 10 4 20])
set(gca,'Fontsize',13,'XTick',0:2:15,'YTick',0:1:20)
xlabel('Ocurrences')
ylabel('Bubble diameter (mm)')

[Sm,Ss] = deal(nan(nfm,1));

% Timeseries
axes('Position',[0.04 0.89 0.945 0.1])
hold on
sur = patch([(1:nfm)';(nfm:-1:1)'],[(Sm-Ss);flipud(Sm+Ss)],'r');
set(sur,'FaceColor',[0.5 0.5 0.5],'EdgeColor','None')
stp = plot(1:nfm,nan(nfm,1),'k','Linewidth',2);
axis([0 nfm 5 15])
set(gca,'Fontsize',13,'YTick',0:2:20)
xlabel('Frames')
ylabel('Mean \phi (mm)')
grid on; box on
legend({'Mean value','Standard deviation'},'Location','NE','Orientation','Horizontal')

% Recording the algorithm
% film = VideoWriter('Bubbles.mp4','MPEG-4');
% film.FrameRate = 24;
% film.Quality = 100;
% open(film);

SS = cell(nfm,1);
for ff = 1 : nfm

    % Loading frame
    df = 0;
    fra = read(vid,[df+ff df+ff]);
    img = fra(lim(3):lim(4),lim(1):lim(2),:);
    fra = rgb2gray(img);

    % Detecting edges
    edg = edge(fra,'Sobel',0.1);
    CHN = [255 0 0];
    for ll = 1 : 3
        IMG = img(:,:,ll);
        IMG(edg) = CHN(ll);
        img(:,:,ll) = IMG;
    end
    set(pfr,'CData',img)
    
    % Filling bubbles holes
    for ii = 1 : ni
        lin = edg(ii,:);
        if sum(lin) > 1
            I = find(lin==1);
            d = diff(I);
            for jj = 1 : length(d)
                if d(jj)<buf
                    lin(I(jj):I(jj)+d(jj)) = 1;
                end
            end
        end
        edg(ii,:) = lin;
    end
    for ii = 1 : nj
        lin = edg(:,ii);
        if sum(lin) > 1
            I = find(lin==1);
            d = diff(I);
            for jj = 1 : length(d)
                if d(jj)<buf
                    lin(I(jj):I(jj)+d(jj)) = 1;
                end
            end
        end
        edg(:,ii) = lin;
    end

    % Clustering bubbles
    I = find(edg==1);
    G = nan(size(edg));
    [II,JJ] = ind2sub([ni,nj],I);
    K = 1;
    G(I(1)) = K;
    for ii = 2 : length(I)
        i = [I(ii)-1,I(ii)+1,I(ii)-ni,I(ii)+ni];
        [A,B] = ind2sub([ni,nj],i);
       
        % If on the edges
        f1 = find(A<1|A>ni);
        f2 = find(B<1|B>nj);
        if ~isempty(f1) || ~isempty(f2)
            i(f1) = [];
            i(f2) = [];
        end
    
        p = edg(i).*G(i);
        if sum(p,'omitnan')>0
            j = find(p>0,1,'first');
            G(I(ii)) = p(j);
        else
            K = K+1;
            G(I(ii)) = K;
        end
    end

    % Reducing number of clusters
    L1 = length(unique(G(~isnan(G))));
    L2 = L1+1;
    while L2 > L1
        H = 1:max(G(:));
        for ii = 1 : length(H)
            F = find(G==H(ii));
            for jj = 1 : length(F)
                i = [F(jj)-1,F(jj)+1,F(jj)-ni,F(jj)+ni];
                [A,B] = ind2sub([ni,nj],i);
               
                % If on the edges
                f = A>=1&A<=ni&B>=1&B<=nj;
                i = i(f);

                D = G(i)-G(F(jj));
                f = find(~isnan(D) & D~=0,1,'first');
                if ~isempty(f)
                    G(G==G(i(f))) = G(F(jj));
                end
            end
        
        end
        L2 = length(unique(G(~isnan(G))));
    end

    Q = unique(G(~isnan(G)));
    C = zeros(L2,2);
    [S,O,L] = deal(zeros(L2,1));
    for ii = 1 : L2
        I = find(G==Q(ii));
        [A,B] = ind2sub([ni,nj],I);
        C(ii,:) = [mean(B),mean(A)];
        L(ii) = length(I); % Number of pixels
        O(ii) = res*L(ii); % Total area
        S(ii) = 2*sqrt(res*length(I)/pi); % Bubbles diameter
    end

    % Filtering bubbles smaller than 2x2 pixels
    I = L>3;
    O = O(I);
    C = C(I,:);
    S = S(I);
    SS{ff} = S;
    [hh,~] = histcounts(S,hi);
    Sm(ff) = mean(S);
    Ss(ff) = std(S);

    set(plt,'XData',C(:,1),'YData',C(:,2),'SizeData',4*S)
    set(brh,'YData',hh)
    set(sur,'XData',[(1:ff)';(ff:-1:1)'],'YData',...
        [Sm(1:ff)-Ss(1:ff);flipud(Sm(1:ff)+Ss(1:ff))])
    set(stp,'XData',1:nfm,'YData',Sm)
    drawnow
    % writeVideo(film,getframe(gcf));
end
% close(film)

%% Some statistics
B = cell2mat(SS);
Bm = mean(B);
Bd = median(B);
AA = 1e-6*res*numel(fra);
Ds = zeros(size(SS));
for ii = 1 : length(SS)
    Ds(ii) = length(SS{ii});
end
Ds = Ds./AA;
ih = 1:2:61;
Dm = mean(Ds);
Dd = median(Ds);

figure('Position',[50 50 1280 800],'Color','w')
axes('Position',[0.075 0.59 0.91 0.39])
hold on
[h1,~] = histcounts(B,hi);
brh1 = bar(hi(1:end-1)+0.5*diff(hi(1:2)),h1);
grid on
box on
xlim([4 22]);
set(gca,'Fontsize',18,'XTick',0:1:22,'Layer','Top')
xlabel('Bubble diameter (mm)')
ylabel('Occurrences')

axes('Position',[0.075 0.1 0.91 0.39])
hold on
[h2,~] = histcounts(Ds,ih);
brh2 = bar(ih(1:end-1)+0.5*diff(ih(1:2)),h2);
grid on
box on
xlim([3 39]);
set(gca,'Fontsize',18,'XTick',0:2:100,'Layer','Top')
xlabel('Bubbles/m^2')
ylabel('Occurrences')



















