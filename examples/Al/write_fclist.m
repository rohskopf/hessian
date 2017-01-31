%%% Use the neighlist file and Hessian to write a FC list

h = dlmread('Al222.hessian');
nl = dlmread('NEIGHLIST');
N = 32;
fh_ii = fopen('IILIST', 'w');
fh_ij = fopen('IJLIST', 'w');

% Now need to make an array fclist that has the same order of itag jtag as neighlist
[length_h width_h] = size(h);
[length_nl width_nl] = size(nl);

%fputs(fh_ii, [num2str(length_nl) "\n"]);
fputs(fh_ii, "ii\n");
fputs(fh_ij, "ij\n");

iilist = [];
for a = 1:N
  
  % i-i interaction
  indices = find(h(:,1) == a & h(:,3) == a);
  fcs = h(indices, :);
  iilist = [iilist; fcs];

end

ijlist = [];
for a=1:length_nl

  l = nl(a,:);
  itag = l(1);
  jtag = l(2);

  % i-j interaction
  indices = find(h(:,1) == itag & h(:,3) == jtag);
  fcs = h(indices, :);
  vals = fcs(:,5);

  % Apply multiplicity operations
  %if (any(abs(vals) == abs(-0.0002764201604633054) ) )
  %  disp("yes\n");
  %endif


  ind = abs(vals) == abs(-4.460156772762287e-04);
  vals(ind) = vals(ind)./2;
  ind = abs(vals) == abs(-1.244234793234594e-03);
  vals(ind) = vals(ind)./2;
  ind = abs(vals) == abs(7.225148813893641e-05);
  vals(ind) = vals(ind)./4;
  ind = abs(vals) == abs(7.770240557741474e-03);
  vals(ind) = vals(ind)./4;
  ind = abs(vals) == abs(-1.488103441587289e-05);
  vals(ind) = vals(ind)./2;
  ind = abs(vals) == abs(-1.178995544337022e-04);
  vals(ind) = vals(ind)./2;
  ind = abs(vals) == abs(-7.231376457584979e-04);
  vals(ind) = vals(ind)./2;
  ind = abs(vals) == abs(-1.290522497922984e-03);
  vals(ind) = vals(ind)./2;
  %ind = abs(vals) == abs(-4.298283445800675e-04);
  %vals(ind) = vals(ind)./8;

  %ind = abs(vals) == abs(-6.407223278066353e-02);
  %vals(ind) = vals(ind)./4;
  %ind = abs(vals) == abs(-4.327567893296926e-02);
  %vals(ind) = vals(ind)./4;





  %ind = abs(vals) == abs(-2.764201604633054e-04);
  %vals(ind) = vals(ind)./2;
  %ind = abs(vals) == abs(-0.001074639276421671);
  %vals(ind) = vals(ind)./2;
  %ind = abs(vals) == abs(1.547144823970501e-04);
  %vals(ind) = vals(ind)./2;
  %ind = abs(vals) == abs(1.178995544337022e-04);
  %vals(ind) = vals(ind)./2;
  %ind = abs(vals) == abs(-5.535421289455759e-04);
  %vals(ind) = vals(ind)./2;
  %ind = abs(vals) == abs(1.178995544337022e-04);
  %vals(ind) = vals(ind)./2;
  %ind = abs(vals) == abs(1.290522497922981e-03);
  %vals(ind) = vals(ind)./2;
  %ind = abs(vals) == abs(1.290522497922984e-03);
  %vals(ind) = vals(ind)./2;

  %ind = vals == -4.298283445800675e-04;
  %vals(ind) = vals(ind)./8;

  %ind = vals == -9.337491102564884e-04;
  %vals(ind) = vals(ind)./2;
  %ind = vals == -1.731968226214851e-03;
  %vals(ind) = vals(ind)./2;


  fcs(:,5) = vals;
  ijlist = [ijlist; fcs];

end


% Write the fclist to a file
dlmwrite(fh_ii, iilist, ' ');
dlmwrite(fh_ij, ijlist, ' ');

fclose(fh_ii);
fclose(fh_ij);

