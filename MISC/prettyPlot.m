%make a pretty plot
%                              \
%                               \
%                                \\
%                                 \\
%                                  >\/7
%                              _.-(6'  \
%                             (=___._/` \
%                                  )  \ |
%                                 /   / |
%                                /    > /   
%                               j    < _\
%                           _.-' :      ``.
%                           \ r=._\        `.
%                          <`\\_  \         .`-.
%                           \ r-7  `-. ._  ' .  `\
%                            \`,      `-.`7  7)   )
%                             \/         \|  \'  / `-._
%                                        ||    .'
%                                         \\  (
%                                          >\  >
%                                      ,.-' >.'
%                                     <.'_.''
%                                       <'

%%
%plot whatever you little heart desires here
figure(1)

subplot(2,1,1)
plot(f,20*log10(abs(h)),'LineWidth',2)
grid on;

%%
%make that plot look like a freaking movie star! 
ax = gca;
set(gca,'FontSize',15);
set(gca,'LineWidth',1);
ax.YLim = [-100 20];

%%
%change these labels to fit your needs
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
suptitle('Magnitude and Phase Response of Low Pass Filter Used to Remove Artifact')

%%
%if it happens to be a subplot just do that shiz again
subplot(2,1,2)
plot(f,unwrap(angle(h))*(180/pi),'LineWidth',2)
grid on;
%%
%make that plot look like a freaking movie star! 
set(gca,'FontSize',15);
set(gca,'LineWidth',1);
ax.YLim = [-100 20];

%%
%change these labels to fit your needs
xlabel('Frequency (Hz)')
ylabel('Phase (degrees)')


%Other useful commands

%make plot size of screen

% set(gca,'FontSize',15);
% set(gca,'LineWidth',2);
% set(gcf,'Position',get(0,'Screensize')); 
% %make background white
% 
% set(gcf,'Color','w');


%now all other plots will be so totally jelly of your plot... so jelly they
%will be jam possibly even preservatives. 
% 
%             _____
%            <_____>
%           /  ___  \
%           |-'   `-|
%           |  JAM  |
%           |`-...-'|
%            `=====' 

