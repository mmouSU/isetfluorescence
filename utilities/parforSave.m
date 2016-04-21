function parforSave( fName, varargin )

    for i=1:numel(varargin)
        eval([inputname(i+1) '=varargin{i};']);
    end
    
    save('-mat',fName,inputname(2));
    for i=2:numel(varargin)
        save('-mat',fName,inputname(i+1),'-append');
    end

end

