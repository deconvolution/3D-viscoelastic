function z=vn(i,j,k,l)
% transform c_{i,j,k,l} to C_{i,j}
z=zeros(1,2);
switch i
    case 1
        switch j
            case 1
                z(1)=1;
            case 2
                z(1)=6;
            case 3
                z(1)=5;
        end
    case 2
        switch j
            case 1
                z(1)=6;
            case 2
                z(1)=2;
            case 3
                z(1)=4;
        end
    case 3
        switch j
            case 1
                z(1)=5;
            case 2
                z(1)=4;
            case 3
                z(1)=3;
        end
end

switch k
    case 1
        switch l
            case 1
                z(2)=1;
            case 2
                z(2)=6;
            case 3
                z(2)=5;
        end
    case 2
        switch l
            case 1
                z(2)=6;
            case 2
                z(2)=2;
            case 3
                z(2)=4;
        end
    case 3
        switch l
            case 1
                z(2)=5;
            case 2
                z(2)=4;
            case 3
                z(2)=3;
        end
end

if z(1)-z(2)>0
    z=flip(z,2);
end
end