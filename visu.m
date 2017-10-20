
%---------------------------------------------------
%---------------------------------------------------
%
% Affichage de la sequence d'images. 
%
% le parametre SEQUENCE designe le repertoire ou 
% se trouve la sequence a visualiser
%
% Elise Arnaud
%
%---------------------------------------------------
%---------------------------------------------------


%--------------------------------------------
% definit la sequence d'images de travail
%--------------------------------------------

% type d'image
EXTENSION = 'png';

% repertoire ou est la sequence a traiter
SEQUENCE = '../seq3/';
%SEQUENCE = '../data/';

% numero de la premiere image a traiter
START = 1;



% load le nom des images de la sequence

filenames = dir([SEQUENCE '*.' EXTENSION]);
filenames = sort({filenames.name});
T = length(filenames);

figure;
set(gcf, 'DoubleBuffer', 'on');

% load la premiere image et l'affiche

tt = START;

im = imread([SEQUENCE filenames{tt}]);
imagesc(im);


while tt <= T-1;
    % load l'image suivante
    
    tt = tt+1;    
    
    im = imread([SEQUENCE filenames{tt}]);
    if isempty(im) 
        break;
    end
    
    % affiche l'image et le num de 
    imagesc(im);
    text(10,10,num2str(tt),'Color', 'r');
    disp(['---------- ' filenames{tt} ' ---------------' ]);
    pause;  
  
    
end


close all


