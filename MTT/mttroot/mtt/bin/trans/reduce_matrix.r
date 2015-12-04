     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
     %%%%% Model Transformation Tools %%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reduce function to write matrices in reduce form
% P.J.Gawthrop  January 22 1998
% Copyright (c) P.J.Gawthrop, 1998


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Version control history
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % $Id: reduce_matrix.r,v 1.2 2000/08/24 09:09:22 peterg Exp $
% % $Log: reduce_matrix.r,v $
% % Revision 1.2  2000/08/24 09:09:22  peterg
% % *** empty log message ***
% %
% % Revision 1.1  1998/01/22 10:39:01  peterg
% % Initial revision
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




PROCEDURE Reduce_Matrix;
BEGIN
  IF MTT_Matrix_n>0 THEN
    IF MTT_Matrix_m>0 THEN
    BEGIN
     write "clear ", MTT_Matrix_name, ";";
     write "MATRIX ", MTT_Matrix_name, "(", MTT_Matrix_n, ",",
                                            MTT_Matrix_m, ");";
    
    FOR i := 1:MTT_Matrix_n DO
      FOR j := 1:MTT_Matrix_m DO IF MTT_Matrix(i,j) NEQ 0 THEN
        write MTT_Matrix_name, "(", i, ",", j, ") := ",
                                            MTT_Matrix(i,j), ";";
    END;
END;

END;;