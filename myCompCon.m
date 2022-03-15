% connected component separation of non-zero elements of an image
% parameter: typcon : type of connectedness (4 ou 8)
%            background : number of layers (backgroubd < min(image) => no background)
 

function ImOut = myCompCon(ImIn,typcon,fond)
[nblin,nbcol] = size(ImIn);
ImIn1=ones(nblin+2,nbcol+2)*fond;
ImIn1(2:nblin+1,2:nbcol+1) = ImIn;
ImOut1 = ones(nblin+2,nbcol+2);
nbclasses= 1;
ptrtab=[1:nblin*nbcol];
con8 = (typcon == 8);
for lin = 2:nblin+1
   for col = 2:nbcol+1
      k = ImIn1(lin,col);
      if (k == fond)
         ImOut1(lin,col) = 1;
      else
         val(1) = (ImIn1(lin-1,col-1)==k)* con8;
         val(2) = (ImIn1(lin-1,col)==k);
         val(3) = (ImIn1(lin-1,col+1)==k) * con8;
         val(4) = (ImIn1(lin,col-1)==k);
         lab(1) = ImOut1(lin-1,col-1);
         lab(2) = ImOut1(lin-1,col) ;
         lab(3) = ImOut1(lin-1,col+1);
         lab(4) = ImOut1(lin,col-1);
         lab = lab .* val;
         switch sum(val)
         case 0, 
            nbclasses = nbclasses + 1; ImOut1(lin,col) = nbclasses;
         case 1,
            ImOut1(lin,col) = max(lab);
         otherwise, % handling fusionning array
            lmax = lab(1);
            for i=2:4
               l1 = lmax; l2=lab(i);
               if (min(l1,l2) <= 1)
                  lmax=max(l1,l2);
               else
                  lmax1=l1;
                  while (ptrtab(lmax1)~=lmax1) lmax1=ptrtab(lmax1); end;
                  lmax2=l2;
                  while (ptrtab(lmax2)~=lmax2) lmax2=ptrtab(lmax2); end;
                  lmax = max(lmax1,lmax2);
                  ptrtab(min(lmax1,lmax2)) = lmax;
               end;
            end; %
            ImOut1(lin,col) = lmax;
         end; 
      end; 
   end; 
end; 
% update of ptrtab array

ptrtab1=[1:nbclasses];
nbclasse1 = 0;
for i = 1:nbclasses
   if ptrtab(i)==i
      nbclasse1 = nbclasse1+1;
      ptrtab1(i) = nbclasse1;
   end;
end;
   
for i=nbclasses:-1:1
   ptrtab(i)=ptrtab(ptrtab(i));
end;

for i = 1:nbclasses
   ptrtab(i)= ptrtab1(ptrtab(i));
end;


% update of labels image
ImOut = ImOut1(2:nblin+1,2:nbcol+1);
ImOut = ptrtab(ImOut)-1;

       
         
         
            


