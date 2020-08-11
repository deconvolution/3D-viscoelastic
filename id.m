function Id=id(nx,ny,nz,eqi,i,j,k,tar,i2,j2,k2)
Id=zeros(1,2);
Id(1)=(k-1)*nx*ny*3+((j-1)*nx+i-1)*3+1+eqi-1;
Id(2)=(k2-1)*nx*ny*3+((j2-1)*nx+i2-1)*3+1+tar-1;
end