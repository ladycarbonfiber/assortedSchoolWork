function [] = fracc()
yn = 550;
xn = 963;
xs = linspace(-2.5,1,xn);
ys = linspace(-1,1,yn);
clf;
for(i = 1:xn)
    for(j=1:yn)
        fraccPlot(xs(i), ys(j));
    end
end

end



function [] = fraccPlot( xp,yp )
%FRACC Summary of this function goes here
%   Detailed explanation goes here
x = 0;
y = 0;
hold on;
it = 0;
max = 100;
while((x^2 + y^ 2 < 4) && (it < max))
    xtemp = x^2 - y^2 + xp;
    y = 2 * x * y + yp;
    x = xtemp;
    it = it + 1;
end
color = 'w';
if((x^2 + y^ 2 < 4))
    color = 'k';
else
     if(it < 6)
         color = 'w';
     elseif(it < 14)
            color = 'y';
     elseif(it < 26)
           color = 'r';
     elseif( it < 45)
            color = 'g';
     elseif(it < 100)
            color = 'b';
        
    end
end
plot(xp,yp,'color',color);
end

