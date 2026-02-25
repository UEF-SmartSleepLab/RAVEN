function c = listboxToCell(x)
    % Converts any listbox Items or Value to cell array of char
    if isempty(x)
        c = {};
    elseif ischar(x)
        c = {x};
    elseif isstring(x)
        c = cellstr(x);
    elseif iscell(x)
        c = x;
    else
        error('Unsupported listbox type');
    end
end