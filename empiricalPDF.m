function pdf_value = empiricalPDF(counts,edges,value)

    % Find the bin containing the value
    bin_index = find(value >= edges(1:end-1) & value < edges(2:end));

    % Calculate the empirical PDF value
    if ~isempty(bin_index)
        pdf_value = counts(bin_index);
    else
        pdf_value = 0;  % Value is outside the range of realizations
    end
end