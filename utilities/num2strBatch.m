function out = num2strBatch(num)
  out = cell(1,length(num));
  for i = 1:length(num)
    out{i} = num2str(num(i));
  end