function [deltaPs] = dP(data, numFind)
    data = data(:, 3:13);
    rowsnumFind = find(data.Var3 == numFind);
    data = table2array(data);
    averages = [];
    for i = 1:length(rowsnumFind)
     rowIndex = rowsnumFind(i);
     if rowIndex <= size(data, 1)
     RowAverage = mean(data(rowIndex, 2:end), 'omitnan');
     averages = [averages; RowAverage];
     end
     end
    deltaPs = averages;
end