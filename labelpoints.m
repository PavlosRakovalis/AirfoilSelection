function [h, ext] = labelpoints (xpos, ypos, labels, varargin)
   
if isnumeric(labels) == 1
    labels = num2cell(labels);
end


if ischar(labels)
    labels = cellstr(labels);
end


labels_size = size(labels);
if labels_size(1) >1 && labels_size(2) >1
  
    numericIdx = cellfun(@isnumeric, labels);
    labels(numericIdx) = strsplit(num2str([labels{numericIdx}]));
   
    tempLabs = cell(1, size(labels,1));
    for r = 1:size(labels,1)
        tempLabs{r} = strjoin(labels(r,:), ' ');
    end
    labels = tempLabs;
elseif length(labels_size) > 2
    error('''LABELS'' may be one or two dimensional.')
end


if length(xpos)==1 && length(ypos)>1
    xpos = repmat(xpos, size(ypos));
elseif length(ypos)==1 && length(xpos)>1
    ypos = repmat(ypos, size(xpos));
end


if length(labels)==1 && length(xpos) > 1
    labels = repmat(labels, [1, length(xpos)]);
end


if iscolumn(xpos);      xpos = xpos';       end
if iscolumn(ypos);      ypos = ypos';       end
if iscolumn(labels);    labels = labels';   end


if isequal(length(xpos), length(ypos), length(labels)) == 0 && sum(strcmp('stacked', varargin))==0
    error('xpos, ypos, and labels must all be the same length unless using one input for labels.')
end


if isa(xpos, 'double')
    xinf = find(xpos==inf);
    yinf = find(ypos==inf);
    findinf = [xinf yinf];
    if ~isempty(findinf)
        xpos(findinf)=[];
        ypos(findinf)=[];
        labels(findinf) = [];
    end
end


validPositions = {'N' 'NE' 'E' 'SE' 'S' 'SW' 'W' 'NW' 'center' 'C'};
checkPosition = @(x) any(validatestring(x, validPositions));
checkCoordinate = @(x) (isnumeric(x) | isa(x, 'datetime'));

p = inputParser;
p.FunctionName = mfilename;
addRequired(p, 'xpos', checkCoordinate);
addRequired(p, 'ypos', checkCoordinate);
addRequired(p, 'labels');
addOptional(p, 'position', 'NW', checkPosition);
addOptional(p, 'buffer', 0, @isnumeric);
addOptional(p, 'adjust_axes', 0, @isnumeric);

addParameter(p, 'outliers_SD', 3, @isnumeric);
addParameter(p, 'outliers_Q', 1.5, @isnumeric);
addParameter(p, 'outliers_N', 1, @isnumeric);
addParameter(p, 'outliers_lim', [0,0;0,0]);
addParameter(p, 'outliers_lin', {'SD',1});
addParameter(p, 'outliers_rand', 0.5, @isnumeric);

addParameter(p, 'stacked', 'down');

addParameter(p, 'axHand', 0, @ishandle);
addParameter(p, 'FontSize', 10, @isnumeric);
addParameter(p, 'FontWeight', 'normal');
addParameter(p, 'Color', 'k');
addParameter(p, 'BackgroundColor', 'none');
addParameter(p, 'rotation', 0, @isnumeric);
addParameter(p, 'interpreter', 'none');
parse(p, xpos, ypos, labels, varargin{:})


axHand = p.Results.axHand;
if axHand == 0
    axHand = gca;
end


buf = p.Results.buffer;
xscl = get(axHand, 'xscale');
yscl = get(axHand, 'yscale');
if (~strcmp(xscl, 'linear') || ~strcmp(yscl, 'linear')) && buf~=0
    warning('Buffer size changed to 0 due to non linear axis scales.')
    buf = 0;
end


if isa(xlim(axHand), 'datetime') && ~isa(ylim(axHand), 'datetime')
    bufferUnits = 'x_datetime';
elseif isa(ylim(axHand), 'datetime') && ~isa(xlim(axHand), 'datetime')
    bufferUnits = 'y_datetime';
elseif isa(xlim(axHand), 'datetime') && isa(ylim(axHand), 'datetime')
    bufferUnits = 'xy_datetime';
else
    bufferUnits = 'other';
end


multiColor = false; 
labelColors = p.Results.Color;
nColors = size(labelColors,1);
if isnumeric(labelColors) && nColors>1
    multiColor = true;
    if nColors < length(labels)
        labelColors = repmat(labelColors, ceil(length(labels)/nColors), 1);
    end
end


[va, ha, u1, u2] = get_compass(upper(p.Results.position), buf, bufferUnits);
    function [va, ha, u1, u2] = get_compass(compass_str, buffer, bufferUnits)
        
        switch bufferUnits
            case 'normalize'
                a = [0 1 0 1]/10;
            case 'x_datetime'
                a = [0 0, ylim(axHand)/10];
                if buf~=0 && ismember(compass_str, {'E', 'W', 'NE', 'NW', 'SE', 'SW'})
                    warning('X axis is in datetime units so label buffer must be 0 for East/West orientations.')
                end
            case 'y_datetime'
                a = [0 0, xlim(axHand)/10];
                if buf~=0 && ismember(compass_str, {'N', 'S', 'NE', 'NW', 'SE', 'SW'})
                    warning('Y axis is in datetime units so label buffer must be 0 for North/South orientations.')
                end
            case 'xy_datetime'
                a = [0 0 0 0];
                if buf~=0
                    warning('X & Y axes are in datetime units so label buffer must be 0.')
                end
            otherwise
                a = axis(axHand)/10;
        end
        
        
        u1 = 0; u2 = 0;
        
        switch upper(compass_str)
            case 'E',       va = 'middle'; ha = 'left';         u1 = a(2)-a(1);
            case 'W',       va = 'middle'; ha = 'right';        u1 = (a(2)-a(1))*-1;
            case 'N',       va = 'bottom'; ha = 'center';       u2 = a(4)-a(3);
            case 'S',       va = 'top';    ha = 'center';       u2 = (a(4)-a(3))*-1;
            case 'NE',      va = 'bottom'; ha = 'left';         u1 = (a(2)-a(1))/2;     u2 = (a(4)-a(3))/2;
            case 'NW',      va = 'bottom'; ha = 'right';        u1 = (a(2)-a(1))*-0.5;  u2 = (a(4)-a(3))/2;
            case 'SE',      va = 'top';    ha = 'left';         u1 = (a(2)-a(1))/2;     u2 = (a(4)-a(3))*-0.5;
            case 'SW',      va = 'top';    ha = 'right';        u1 = (a(2)-a(1))*-0.5;  u2 = (a(4)-a(3))*-0.5;
            case {'CENTER', 'C'},  va = 'middle'; ha = 'center';
        end
        
        
        u1 = u1*buffer;   
        u2 = u2*buffer;   
    end

if p.Results.rotation~=0
    factor = 1.2; 
    u1 = u1*factor;
    u2 = u2*factor;
    
    va = 'middle'; ha = 'center';
end


if sum(strcmp('stacked', varargin))==1

    if length(xpos)>1 || length(ypos)>1
        warning('Only the first xpos and ypos will be used to initiate stacked text');
    end
   
    if ~isempty(cell2mat(regexp(varargin(cellfun(@ischar, varargin)), 'outlier')))
        warning('Cannot use stacked and outlier input parameters at the same time.  Outliers input will be ignored');
        tempvarargin = varargin;
        tempvarargin(cellfun(@isnumeric, tempvarargin)) = {'temp_replace'};    
        varargIDX = find(~cellfun(@isempty,regexp(tempvarargin, 'outlier')));
        varargin([varargIDX, varargIDX+1]) = [];    
    end
    
   
    spacing = 1;
    stacktype = lower(p.Results.stacked);
    
 
    if ~isempty(strfind(stacktype, '_'))
        spacing = str2double(stacktype(strfind(stacktype, '_')+1:end));
        stacktype = stacktype(1:strfind(stacktype, '_')-1);
    end
    
 
    if ~any(strcmp(stacktype, {'up', 'down', 'left', 'right'}))
        error('Text stacking options are:  up, down, left, or right (not case sensitive)');
    end
    

    nlabs = length(labels); 
    if nlabs>1            
        xpos(min(2, nlabs):nlabs) = nan;
        ypos(min(2, nlabs):nlabs) = nan;
    end
    

    for s = 2:nlabs
        
     
        labhand = text(xpos(s-1)+u1 , ypos(s-1)+u2, labels(s-1), 'VerticalAlignment',va, 'HorizontalAlignment',ha, 'FontSize', p.Results.FontSize, 'Parent', axHand);
        label_extnt_norm = get(labhand, 'extent');
        delete(labhand)
        
     
        switch stacktype
            case 'down'
                xpos(s) = xpos(s-1);
                ypos(s) = ypos(s-1) - label_extnt_norm(4)*spacing;
            case 'right'
                ypos(s) = ypos(s-1);
                xpos(s) = xpos(s-1) + label_extnt_norm(3)*spacing;
            case 'up'
                xpos(s) = xpos(s-1);
                ypos(s) = ypos(s-1) + label_extnt_norm(4)*spacing;
            case 'left'
                ypos(s) = ypos(s-1);
                xpos(s) = xpos(s-1) - label_extnt_norm(3)*spacing;
        end 
    end 
end 


outlierNames = {'outliers_SD', 'outliers_Q', 'outliers_N', 'outliers_lim', 'outliers_lin', 'outliers_rand'}; 
outlier_flag = false;

for c = 1:length(outlierNames)
    if sum(strcmp(varargin, outlierNames{c}))>0
        outlier_flag = true;
        outliertype = varargin{strcmp(varargin, outlierNames{c})};  
    end
end

if outlier_flag
  
    pairIDX = ~isnan(xpos) & ~isnan(ypos);
    switch outliertype
        
        case 'outliers_SD'
            SDs = p.Results.outliers_SD(1);
           
            if length(p.Results.outliers_SD) > 1
                Xcnt = p.Results.outliers_SD(2);
                Ycnt = p.Results.outliers_SD(3);
            else 
                Xcnt = nanmean(xpos);
                Ycnt = nanmean(ypos);
            end
            outlier_idx = logical(abs(xpos - Xcnt) > SDs*nanstd(xpos)  |  abs(ypos - Ycnt) > SDs*nanstd(ypos)); 
            
        case 'outliers_Q'
            xbounds = [prctile(xpos,25) - p.Results.outliers_Q * iqr(xpos) , prctile(xpos, 75) + p.Results.outliers_Q * iqr(xpos)];   
            ybounds = [prctile(ypos,25) - p.Results.outliers_Q * iqr(ypos) , prctile(ypos, 75) + p.Results.outliers_Q * iqr(ypos)];   
            outlier_idx = logical(ypos<ybounds(1) | ypos>ybounds(2) |  xpos<xbounds(1) | xpos>xbounds(2));
            
        case 'outliers_lim'
           
            limvars = p.Results.outliers_lim; 
            if iscell(limvars)
                lims = limvars{1};
                if size(limvars,1) == 1 
                    qualifier = 'or';
                else
                    qualifier = lower(limvars{2});
                end
            else
                lims = limvars;
                qualifier = 'or';
            end
            
            if size(lims,1) == 1
                lims = [lims;lims];
            end
            
            x_outliers = xpos<lims(1,1) | xpos>lims(1,2);  
            y_outliers = ypos<lims(2,1) | ypos>lims(2,2); 
            switch qualifier
                case 'or'
                    outlier_idx = x_outliers | y_outliers;
                case 'and'
                    outlier_idx = x_outliers & y_outliers;
                case 'invert'
                    x_outliers = xpos>lims(1,1) & xpos<lims(1,2);
                    y_outliers = ypos>lims(2,1) & ypos<lims(2,2);
                    outlier_idx = x_outliers & y_outliers;
                otherwise
                    error('The inputs you entered for Outliers_lim wasn''t recognized.')
            end
            
        case 'outliers_N'
            Npts = p.Results.outliers_N(1);
           
            if length(p.Results.outliers_N) > 1
                Xcnt = p.Results.outliers_N(2);
                Ycnt = p.Results.outliers_N(3);
            else 
                Xcnt = nanmean(xpos);
                Ycnt = nanmean(ypos);
            end
            if p.Results.outliers_N<1
                N = ceil(length(xpos(pairIDX)) * Npts);
            else
                N = min(Npts, length(xpos(pairIDX)));     
            end
            meanpoint = repmat([Xcnt Ycnt], [length(xpos),1]);
            paired = horzcat(xpos', ypos');
            distances = (((meanpoint(:,1)-paired(:,1)).^2)  +  ((meanpoint(:,2)-paired(:,2)).^2)).^(1/2);      
            [sorted, idx] = sort(distances, 'descend');
            idx = idx(~isnan(sorted));
            outlier_idx = false(1,length(xpos));
            outlier_idx(idx(1:N))=1;
            
        case 'outliers_lin'
            
         
            linvars = p.Results.outliers_lin;
            if isnumeric(linvars{1}) && ~isempty(linvars{1}) 
                slope = linvars{1};
                yint = linvars{2};
                outtype = upper(linvars{3}); 
                outthresh = linvars{4};      
            else
               
                if isempty(linvars{1}); linvars = linvars([3,4]); end
                slope = nansum((xpos-nanmean(xpos)).*(ypos-nanmean(ypos))) / nansum((xpos-nanmean(xpos)).^2);
                yint = nanmean(ypos) - (slope * nanmean(xpos));                    
                outtype = upper(linvars{1});
                outthresh = linvars{2};
            end
           
            Yest = slope*xpos + yint;
            resid = (Yest - ypos).^2;
          
            [sorted, idx] = sort(resid, 'descend');
            idx = idx(~isnan(sorted)); 
         
            if strcmp(outtype, 'SD')
                outlier_idx = idx(1:sum(nanstd(sorted)*outthresh <= sorted));    
            elseif strcmp(outtype, 'N')
                if outthresh<1
                    N = ceil(length(idx) * outthresh);
                else
                    N = min(outthresh, length(idx));        
                end
                outlier_idx = idx(1:N);
            end
            
        case 'outliers_rand'
           
            Npts = p.Results.outliers_rand(1);
            if p.Results.outliers_rand<1
                N = ceil(length(xpos(pairIDX)) * Npts);
            else
                N = min(Npts, length(xpos(pairIDX)));      
            end
            
            outlier_idx = randsample(length(xpos(pairIDX)),N);
            
            
    end 
    
    xpos = xpos(outlier_idx);
    ypos = ypos(outlier_idx);
    labels = labels(outlier_idx);
    
  
end 




if multiColor
    hand = zeros(size(labels));
    for k = 1:length(labels)
        hand(k) = text(xpos(k)+u1 , ypos(k)+u2, labels{k}, 'VerticalAlignment',va, 'HorizontalAlignment',ha, 'FontSize', p.Results.FontSize, 'color', labelColors(k,:), 'FontWeight', p.Results.FontWeight, 'Parent', axHand);
    end
else
    
    hand = text(xpos+u1 , ypos+u2, labels, 'VerticalAlignment',va, 'HorizontalAlignment',ha, 'FontSize', p.Results.FontSize, 'color', labelColors, 'BackgroundColor', p.Results.BackgroundColor, 'FontWeight', p.Results.FontWeight, 'Parent', axHand);
end
extnt = get(hand, 'extent');


if p.Results.rotation~=0           
    xl = xlim;      yl = ylim;                          
    curr_extent = get(hand, 'extent');                     
    if iscell(curr_extent); curr_extent = cell2mat(curr_extent); end
    hold on
    curr_position = [curr_extent(:,1)+(curr_extent(:,3)/2),curr_extent(:,2)+(curr_extent(:,4)/2)];          
    set(hand, 'rotation', p.Results.rotation, 'VerticalAlignment','middle', 'HorizontalAlignment','center');  	
    for i = 1:length(hand)                                
        set(hand(i), 'position', curr_position(i,:))
    end
    set(axHand, 'xlim', xl); set(axHand, 'ylim', yl);         
end


if p.Results.adjust_axes == 1   &&   ~isempty(hand)
    x_adj = sign(u1+0.0000001);                
    y_adj = sign(u2+0.0000001);              
    
    labelextent = get(hand, 'extent');
    if isequal(class(labelextent),'cell')
        labelextent = cat(1, labelextent{:});
    end
    xl = xlim;      yl = ylim;
    lablimX = [min(labelextent(:,1)), max(labelextent(:,1)+(labelextent(:,3).*x_adj))] +u1;
    lablimY = [min(labelextent(:,2)), max(labelextent(:,2)+(labelextent(:,4).*y_adj))] +u2;
    
    xlim([min(min(xl), min(lablimX)), max(max(xl), max(lablimX))])
    ylim([min(min(yl), min(lablimY)), max(max(yl), max(lablimY))])
   
end


set(hand, 'interpreter', p.Results.interpreter)


if nargout > 0
    h   = hand;
    ext = extnt;
end

end


