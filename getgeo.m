%funtion creates txt-file readable with gmsh
function getgeo(c,TRI,x,y,v1,v2,v3,e1,e2,e3,ie1)
%parameter
%c - Coefficient vector of spline
%TRI - ConnectivityList of the vertices
% x,y - cartesian coordinates of vertices
% v1,v2,v3 - Lists of indizes of the vertices in each triangle
% e1,e2,e3 - Lists of indizes of the edges of each triangle
% ie1 - list of indizes of the starting vertex of each edge

%barycentric grid on a triangle
B=[1,0,0;0.5,0.5,0;0.5,0,0.5;0,1,0;0,0.5,0.5;0,0,1];

%number of vertices in given triangulation
nv=length(x);

%number of edges in given triangulation
ne=length(ie1);

%number of triangles in given triangulation
nt=size(TRI,1);

%save for every dot of barycentric grid if already in our list
%degree of spline -> 2 -> number of dots = nt+ne
%inlist=zeros(nv+ne,1);

%all the points in barycentric grid -> first column=index, other columns ->
%cartesian coordinates
%includes also help points for bezier curve 
poi=zeros(nv+2*ne,4);


%indizes of points on the edges of triangulation
%5th coulumn -> 1 = we didn't have this edge in our list, -1 = we have it
%in our list
elist=zeros(ne,5);
elist(:,5)=ones(ne,1);

%indizes of edges of each triangle
surflist=zeros(nt,4);

%TRIANG=triangulation(TRI,[x,y]);

%loop over all triangles -> extract coefficients and find values of spline
for i=1:nt
    %extract indizes of coefficients of the polynomial on current
    %triangle
    ind=finde_ind(i,nv,v1,v2,v3,e1,e2,e3);
    indhelp=[ind(2)+ne;ind(3)+ne;ind(5)+ne];

    %define the coefficients on the current triangle
    c_pol=c(ind);
    
    %include the points in our list for gmsh
    poi(ind,1)=ind;
    poi(indhelp,1)=indhelp;
    
    %vertices of the current triangle
    p1=[x(TRI(i,1)),y(TRI(i,1))];
    p2=[x(TRI(i,2)),y(TRI(i,2))];
    p3=[x(TRI(i,3)),y(TRI(i,3))];
    
    %compute cartesian coordinates of the barycentric grid points on cur.
    %triangle 
    %x1=vertex p1, others lexicographic order
    [x1,x2,x3,x4,x5,x6]=quad_nodes_lexi(p1,p2,p3);
    poi(ind,2:3)=[x1;x2;x3;x4;x5;x6];
    poi(indhelp,2:3)=[x2;x3;x5];
    %z-coordinates by using DeCasteljau algorithm
    poi(ind,4)=DeCasteljau(B(:,1),B(:,2),B(:,3),c_pol');
    poi(indhelp,4)=[-0.5*poi(ind(1),4)+2*poi(ind(2),4)-0.5*poi(ind(4),4);...
                    -0.5*poi(ind(1),4)+2*poi(ind(3),4)-0.5*poi(ind(6),4);...
                    -0.5*poi(ind(4),4)+2*poi(ind(5),4)-0.5*poi(ind(6),4)];

    %indizes of edges on current triangle
    e_ind=[e1(i),e2(i),e3(i)];
            
    %for getting CurveLoop on current triangle -> edge indizes with
    %orientation ( if 5th column of elist = -1-> reverse orientation)
    surflist(i,:)=[i,elist(e1(i),5)*e1(i),elist(e2(i),5)*e2(i),elist(e3(i),5)*e3(i)];
    
    %if edges are not in our list -> include
    if elist(e1(i),5) == 1
    elist(e1(i),:)=[e_ind(1),ind(1),ind(2)+ne,ind(4),-elist(e1(i),5)];
    end
    
    if elist(e2(i),5) == 1
    elist(e2(i),:)=[e_ind(2),ind(4),ind(5)+ne,ind(6),-elist(e2(i),5)];
    end
    
    if elist(e3(i),5) == 1
    elist(e3(i),:)=[e_ind(3),ind(6),ind(3)+ne,ind(1),-elist(e3(i),5)];
    end
end


% GMsh geometry file is stored as .txt file
fileID = fopen('example.geo','w');

fprintf(fileID,'// OpenCASCADE suffers with some bugs:\n');
fprintf(fileID,'// https://gitlab.onelab.info/gmsh/gmsh/issues/594\n');
fprintf(fileID,'// https://github.com/tpaviot/oce/issues/716\n');
fprintf(fileID,'//SetFactory("OpenCASCADE");\n');
fprintf(fileID,'\n');
fprintf(fileID,'SetFactory("Built-in");\n');
fprintf(fileID,'\n');

%poi(:, 4) = 0.1*poi(:, 4);

%list of points
for i=1:size(poi,1)
fprintf(fileID,'Point(%d)={%d,%d,%d};\n',poi(i,:));
%fprintf(fileID,'\n');
end

%list of splines
for i = 1:ne
fprintf(fileID,'Bezier(%d) = {%d,%d,%d};\n',elist(i,1:4));
%fprintf(fileID,'Bezier(%d) = {%d,%d,%d};\n',elist(i,1:4));
%fprintf(fileID,'\n');
end

%list of Curve Loops
for i=1:nt
fprintf(fileID,'Curve Loop(%d) = {%d,%d,%d};\n',surflist(i,1:4));
%fprintf(fileID,'\n');
end

%list of surfaces
for i=1:nt
    fprintf(fileID,'Surface(%d)={%d};\n',[i,i]);
    %fprintf(fileID,'\n');
end

%save e.g. brep-File
%fprintf(fileID,'Save "geotest.brep";');

fprintf(fileID,'\n');
fprintf(fileID,'Physical Surface("terrain") = {Surface{:}};\n');
fprintf(fileID,'Compound Surface{Surface{:}};\n');
fprintf(fileID,'\n');
fprintf(fileID,'//Mesh.CharacteristicLengthMin = 50;\n');
fprintf(fileID,'//Mesh.CharacteristicLengthMax = 50;\n');
fprintf(fileID,'\n');
fprintf(fileID,'Mesh 2;\n');
fprintf(fileID,'Save ''example.vtk'';\n');

fclose(fileID);
%type geotest.txt
%uiopen('C:\Users\max-k\Desktop\Bachelorarbeit\splinepak\Computerpraktikum\geotest.txt',1)

end
