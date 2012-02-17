function lineSpec = giveMeLineSpec(ix)
%
%
%

% color
x = rem(ix,7);
switch x
   case 1
      col = 'r';
   case 2
      col = 'g';
   case 3
      col = 'b';
   case 4
      col = 'c';
   case 5
      col = 'm';
   case 6
      col = 'y';
   otherwise
      col = 'k';
end
ix = fix(ix / 8);

% line
x = rem(ix,4);
switch x
   case 0
      lin = '-';
   case 1
      lin = '--';
   case 2
      lin = ':';
   otherwise
      lin = '--';
end
ix = fix(ix/5);

% marker
x = rem(ix,13);
switch x
   case 0
      mar = '';
   case 1
      mar = '+';
   case 2
      mar = 'o';
   case 3
      mar = '*';
   case 4
      mar = '.';
   case 5
      mar = 'x';
   case 6
      mar = 's';
   case 7
      mar = 'd';
   case 8
      mar = '^';
   case 9
      mar = 'v';
   case 10
      mar = '<';
   case 11
      mar = '>';
   case 12
      mar = 'p';
   otherwise
      lin = 'h';
end

%
lineSpec = [lin col mar];