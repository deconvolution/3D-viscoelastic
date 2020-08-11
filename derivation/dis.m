eqd=subs(eq,[... pux
    diff(ux,x,1), ...
    diff(ux,y,1), ...
    diff(ux,z,1), ...
    diff(ux,t,1), ...
    ... p^2 ux
    diff(ux,x,2), ...
    diff(ux,y,2), ...
    diff(ux,z,2), ...
    diff(ux,t,2), ...
    diff(diff(ux,y,1),z,1), ...
    diff(diff(ux,x,1),z,1), ...
    diff(diff(ux,x,1),y,1), ...
    diff(diff(ux,t,1),x,1), ...
    diff(diff(ux,t,1),y,1), ...
    diff(diff(ux,t,1),z,1), ...
    ... p^3 ux
    diff(diff(ux,x,2),t), ...
    diff(diff(ux,y,2),t), ...
    diff(diff(ux,z,2),t), ...
    diff(diff(diff(ux,t,1),y,1),z,1), ...
    diff(diff(diff(ux,t,1),x,1),z,1),...
    diff(diff(diff(ux,t,1),x,1),y,1),...
     ... puy
    diff(uy,x,1), ...
    diff(uy,y,1), ...
    diff(uy,z,1), ...
    diff(uy,t,1), ...
    ... p^2 uy
    diff(uy,x,2), ...
    diff(uy,y,2), ...
    diff(uy,z,2), ...
    diff(uy,t,2), ...
    diff(diff(uy,y,1),z,1), ...
    diff(diff(uy,x,1),z,1), ...
    diff(diff(uy,x,1),y,1), ...
    diff(diff(uy,t,1),x,1), ...
    diff(diff(uy,t,1),y,1), ...
    diff(diff(uy,t,1),z,1), ...
    ... p^3 uy
    diff(diff(uy,x,2),t),...
    diff(diff(uy,y,2),t),...
    diff(diff(uy,z,2),t),...
    diff(diff(diff(uy,t,1),y,1),z,1),...
    diff(diff(diff(uy,t,1),x,1),z,1),...
    diff(diff(diff(uy,t,1),x,1),y,1),...
    ... puz
    diff(uz,x,1), ...
    diff(uz,y,1), ...
    diff(uz,z,1), ...
    diff(uz,t,1), ...
    ... p^2 uz
    diff(uz,x,2), ...
    diff(uz,y,2), ...
    diff(uz,z,2), ...
    diff(uz,t,2), ...
    diff(diff(uz,y,1),z,1), ...
    diff(diff(uz,x,1),z,1), ...
    diff(diff(uz,x,1),y,1), ...
    diff(diff(uz,t,1),x,1), ...
    diff(diff(uz,t,1),y,1), ...
    diff(diff(uz,t,1),z,1), ...
    ... p^3 uz
    diff(diff(uz,x,2),t),...
    diff(diff(uz,y,2),t),...
    diff(diff(uz,z,2),t),...
    diff(diff(diff(uz,t,1),y,1),z,1),...
    diff(diff(diff(uz,t,1),x,1),z,1),...
    diff(diff(diff(uz,t,1),x,1),y,1),...
    ], ...
    ...
    [...pux
    (uxd(4,3,3,3)-uxd(2,3,3,3))/2/dx, ...
    (uxd(3,4,3,3)-uxd(3,2,3,3))/2/dy, ...
    (uxd(3,3,4,3)-uxd(3,3,2,3))/2/dz, ...
    (uxd(3,3,3,4)-uxd(3,3,3,3))/dt, ...
    ...p^2 ux
    (uxd(4,3,3,3)-2*uxd(3,3,3,3)+uxd(2,3,3,3))/dx^2, ...
    (uxd(3,4,3,3)-2*uxd(3,3,3,3)+uxd(3,2,3,3))/dy^2, ...
    (uxd(3,3,4,3)-2*uxd(3,3,3,3)+uxd(3,3,2,3))/dz^2, ...
    (uxd(3,3,3,4)-2*uxd(3,3,3,3)+uxd(3,3,3,2))/dt^2, ...
    (uxd(3,4,4,3)-uxd(3,2,4,3)-uxd(3,4,2,3)+uxd(3,2,2,3))/4/dy/dz,...
    (uxd(4,3,4,3)-uxd(2,3,4,3)-uxd(4,3,2,3)+uxd(2,3,2,3))/4/dx/dz,...
    (uxd(4,4,3,3)-uxd(2,4,3,3)-uxd(4,2,3,3)+uxd(2,2,3,3))/4/dx/dy,...
    (uxd(4,3,3,3)-uxd(2,3,3,3)-uxd(4,3,3,2)+uxd(2,3,3,2))/2/dt/dx,...
    (uxd(3,4,3,3)-uxd(3,2,3,3)-uxd(3,4,3,2)+uxd(3,2,3,2))/2/dt/dy,...
    (uxd(3,3,4,3)-uxd(3,3,2,3)-uxd(3,3,4,2)+uxd(3,3,2,2))/2/dt/dz,...
    ...p^3 ux
    (uxd(4,3,3,3)-2*uxd(3,3,3,3)+uxd(2,3,3,3))/dt/dx^2-(uxd(4,3,3,2)-2*uxd(3,3,3,2)+uxd(2,3,3,2))/dt/dx^2,...
    (uxd(3,4,3,3)-2*uxd(3,3,3,3)+uxd(3,2,3,3))/dt/dy^2-(uxd(3,4,3,2)-2*uxd(3,3,3,2)+uxd(3,2,3,2))/dt/dy^2,...
    (uxd(3,3,4,3)-2*uxd(3,3,3,3)+uxd(3,3,2,3))/dt/dz^2-(uxd(3,3,4,2)-2*uxd(3,3,3,2)+uxd(3,3,2,2))/dt/dz^2,...
    (uxd(3,4,4,3)-uxd(3,2,4,3)-uxd(3,4,2,3)+uxd(3,2,2,3))/4/dt/dy/dz-(uxd(3,4,4,2)-uxd(3,2,4,2)-uxd(3,4,2,2)+uxd(3,2,2,2))/4/dt/dy/dz,...
    (uxd(4,3,4,3)-uxd(2,3,4,3)-uxd(4,3,2,3)+uxd(2,3,2,3))/4/dt/dx/dz-(uxd(4,3,4,2)-uxd(2,3,4,2)-uxd(4,3,2,2)+uxd(2,3,2,2))/4/dt/dx/dz,...
    (uxd(4,4,3,3)-uxd(2,4,3,3)-uxd(4,2,3,3)+uxd(2,2,3,3))/4/dt/dx/dy-(uxd(4,4,3,2)-uxd(2,4,3,2)-uxd(4,2,3,2)+uxd(2,2,3,2))/4/dt/dx/dy,...
     ...puy
    (uyd(4,3,3,3)-uyd(2,3,3,3))/2/dx, ...
    (uyd(3,4,3,3)-uyd(3,2,3,3))/2/dy, ...
    (uyd(3,3,4,3)-uyd(3,3,2,3))/2/dz, ...
    (uyd(3,3,3,4)-uyd(3,3,3,3))/dt, ...
    ...p^2 uy
    (uyd(4,3,3,3)-2*uyd(3,3,3,3)+uyd(2,3,3,3))/dx^2, ...
    (uyd(3,4,3,3)-2*uyd(3,3,3,3)+uyd(3,2,3,3))/dy^2, ...
    (uyd(3,3,4,3)-2*uyd(3,3,3,3)+uyd(3,3,2,3))/dz^2, ...
    (uyd(3,3,3,4)-2*uyd(3,3,3,3)+uyd(3,3,3,2))/dt^2, ...
    (uyd(3,4,4,3)-uyd(3,2,4,3)-uyd(3,4,2,3)+uyd(3,2,2,3))/4/dy/dz,...
    (uyd(4,3,4,3)-uyd(2,3,4,3)-uyd(4,3,2,3)+uyd(2,3,2,3))/4/dx/dz,...
    (uyd(4,4,3,3)-uyd(2,4,3,3)-uyd(4,2,3,3)+uyd(2,2,3,3))/4/dx/dy,...
    (uyd(4,3,3,3)-uyd(2,3,3,3)-uyd(4,3,3,2)+uyd(2,3,3,2))/2/dt/dx,...
    (uyd(3,4,3,3)-uyd(3,2,3,3)-uyd(3,4,3,2)+uyd(3,2,3,2))/2/dt/dy,...
    (uyd(3,3,4,3)-uyd(3,3,2,3)-uyd(3,3,4,2)+uyd(3,3,2,2))/2/dt/dz,...
    ...p^3 uy
    (uyd(4,3,3,3)-2*uyd(3,3,3,3)+uyd(2,3,3,3))/dt/dx^2-(uyd(4,3,3,2)-2*uyd(3,3,3,2)+uyd(2,3,3,2))/dt/dx^2,...
    (uyd(3,4,3,3)-2*uyd(3,3,3,3)+uyd(3,2,3,3))/dt/dy^2-(uyd(3,4,3,2)-2*uyd(3,3,3,2)+uyd(3,2,3,2))/dt/dy^2,...
    (uyd(3,3,4,3)-2*uyd(3,3,3,3)+uyd(3,3,2,3))/dt/dz^2-(uyd(3,3,4,2)-2*uyd(3,3,3,2)+uyd(3,3,2,2))/dt/dz^2,...
    (uyd(3,4,4,3)-uyd(3,2,4,3)-uyd(3,4,2,3)+uyd(3,2,2,3))/4/dt/dy/dz-(uyd(3,4,4,2)-uyd(3,2,4,2)-uyd(3,4,2,2)+uyd(3,2,2,2))/4/dt/dy/dz,...
    (uyd(4,3,4,3)-uyd(2,3,4,3)-uyd(4,3,2,3)+uyd(2,3,2,3))/4/dt/dx/dz-(uyd(4,3,4,2)-uyd(2,3,4,2)-uyd(4,3,2,2)+uyd(2,3,2,2))/4/dt/dx/dz,...
    (uyd(4,4,3,3)-uyd(2,4,3,3)-uyd(4,2,3,3)+uyd(2,2,3,3))/4/dt/dx/dy-(uyd(4,4,3,2)-uyd(2,4,3,2)-uyd(4,2,3,2)+uyd(2,2,3,2))/4/dt/dx/dy,...
    ...puz
    (uzd(4,3,3,3)-uzd(2,3,3,3))/2/dx, ...
    (uzd(3,4,3,3)-uzd(3,2,3,3))/2/dy, ...
    (uzd(3,3,4,3)-uzd(3,3,2,3))/2/dz, ...
    (uzd(3,3,3,4)-uzd(3,3,3,3))/dt, ...
    ...p^2 uz
    (uzd(4,3,3,3)-2*uzd(3,3,3,3)+uzd(2,3,3,3))/dx^2, ...
    (uzd(3,4,3,3)-2*uzd(3,3,3,3)+uzd(3,2,3,3))/dy^2, ...
    (uzd(3,3,4,3)-2*uzd(3,3,3,3)+uzd(3,3,2,3))/dz^2, ...
    (uzd(3,3,3,4)-2*uzd(3,3,3,3)+uzd(3,3,3,2))/dt^2, ...
    (uzd(3,4,4,3)-uzd(3,2,4,3)-uzd(3,4,2,3)+uzd(3,2,2,3))/4/dy/dz,...
    (uzd(4,3,4,3)-uzd(2,3,4,3)-uzd(4,3,2,3)+uzd(2,3,2,3))/4/dx/dz,...
    (uzd(4,4,3,3)-uzd(2,4,3,3)-uzd(4,2,3,3)+uzd(2,2,3,3))/4/dx/dy,...
    (uzd(4,3,3,3)-uzd(2,3,3,3)-uzd(4,3,3,2)+uzd(2,3,3,2))/2/dt/dx,...
    (uzd(3,4,3,3)-uzd(3,2,3,3)-uzd(3,4,3,2)+uzd(3,2,3,2))/2/dt/dy,...
    (uzd(3,3,4,3)-uzd(3,3,2,3)-uzd(3,3,4,2)+uzd(3,3,2,2))/2/dt/dz,...
    ...p^3 uz
    (uzd(4,3,3,3)-2*uzd(3,3,3,3)+uzd(2,3,3,3))/dt/dx^2-(uzd(4,3,3,2)-2*uzd(3,3,3,2)+uzd(2,3,3,2))/dt/dx^2,...
    (uzd(3,4,3,3)-2*uzd(3,3,3,3)+uzd(3,2,3,3))/dt/dy^2-(uzd(3,4,3,2)-2*uzd(3,3,3,2)+uzd(3,2,3,2))/dt/dy^2,...
    (uzd(3,3,4,3)-2*uzd(3,3,3,3)+uzd(3,3,2,3))/dt/dz^2-(uzd(3,3,4,2)-2*uzd(3,3,3,2)+uzd(3,3,2,2))/dt/dz^2,...
    (uzd(3,4,4,3)-uzd(3,2,4,3)-uzd(3,4,2,3)+uzd(3,2,2,3))/4/dt/dy/dz-(uzd(3,4,4,2)-uzd(3,2,4,2)-uzd(3,4,2,2)+uzd(3,2,2,2))/4/dt/dy/dz,...
    (uzd(4,3,4,3)-uzd(2,3,4,3)-uzd(4,3,2,3)+uzd(2,3,2,3))/4/dt/dx/dz-(uzd(4,3,4,2)-uzd(2,3,4,2)-uzd(4,3,2,2)+uzd(2,3,2,2))/4/dt/dx/dz,...
    (uzd(4,4,3,3)-uzd(2,4,3,3)-uzd(4,2,3,3)+uzd(2,2,3,3))/4/dt/dx/dy-(uzd(4,4,3,2)-uzd(2,4,3,2)-uzd(4,2,3,2)+uzd(2,2,3,2))/4/dt/dx/dy,...
    ]);
eqd=subs(eqd,...
    [...pUx
    diff(Ux,x,1),...
    diff(Ux,y,1),...
    diff(Ux,z,1),...
    ...pUy
    diff(Uy,x,1),...
    diff(Uy,y,1),...
    diff(Uy,z,1),...
    ...pUz
    diff(Uz,x,1),...
    diff(Uz,y,1),...
    diff(Uz,z,1)...
    ],...
    [...pUx
    (Uxd(4,3,3)-Uxd(2,3,3))/2/dx,...
    (Uxd(3,4,3)-Uxd(3,2,3))/2/dy,...
    (Uxd(3,3,4)-Uxd(3,3,2))/2/dz,...
    ...pUx
    (Uyd(4,3,3)-Uyd(2,3,3))/2/dx,...
    (Uyd(3,4,3)-Uyd(3,2,3))/2/dy,...
    (Uyd(3,3,4)-Uyd(3,3,2))/2/dz,...
    ...pUx
    (Uzd(4,3,3)-Uzd(2,3,3))/2/dx,...
    (Uzd(3,4,3)-Uzd(3,2,3))/2/dy,...
    (Uzd(3,3,4)-Uzd(3,3,2))/2/dz...
    ]);
eqd=subs(eqd,...
    [ux(x,y,z,t),uy(x,y,z,t),uz(x,y,z,t),...
    Ux(x,y,z,t),Uy(x,y,z,t),Uz(x,y,z,t)],...
    [uxd(3,3,3,3),uyd(3,3,3,3),uzd(3,3,3,3),...
    Uxd(3,3,3),Uyd(3,3,3),Uzd(3,3,3)]);