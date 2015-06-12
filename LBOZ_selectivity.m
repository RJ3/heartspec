% function ksi=LBOZ_selectivity(S,pinvS)
% In: calibration matrix S (JxK)
% J is the number of wavelengths (~1044)
% K is the number of references (eg. 9)
% pseudoinverse of S (KxJ)
% Out: LBOZ selectivity vector ksi (1xK)

% load('refs.mat');
load('refs_chess_ferryl.mat')
load('new_wavelengths.mat')

%% Selectivity of Analytes 
% Make sure refs is JxK vector
S=loadrefs(479:734,1:8);
pinvS=pinv(S);
ksi=1./sqrt(sum(S.^2).*sum((pinvS.^2)')); %#ok<UDIM> % ALL REFS
ksi'

% ksi=1./(norm(pinvS(1,:)).*norm(S(:,1))); % SELECT 1 ref at a time



%% Prediction Calculation

wl_start=479;
wl_end=733;

win=DelOD(wl_start:wl_end,x5);


% refs(:,10)=ones;
S=refs(wl_start:wl_end,1:9);
S_plus=pinv(S);
% S_plus(10,:)=ones;
c=S_plus*win;

prediction=(c(1)*S(:,1))+(c(2)*S(:,2))+(c(3)*S(:,3))+...
        (c(4)*S(:,4))+(c(5)*S(:,5))+(c(6)*S(:,6))+...
        (c(7)*S(:,7))+(c(8)*S(:,8))+(c(9)*S(:,9));

waxis=wavelengths(wl_start:wl_end);    

figure
plot(waxis,prediction,'k');
hold on
plot(waxis,win,'r');
plot(waxis,win-prediction,'g');
legend('Prediction','Window','Residuals')

%% Standard error of the Prediction

SE=std(win)./(ksi.*sqrt(sum(S.^2))); % 2-norm is the Euclidean norm.
% SE=std(win)./(ksi*norm(S(:,1))); % ONE AT A TIME