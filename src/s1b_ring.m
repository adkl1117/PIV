root = piv_path('data');
C = double(imread([root 'Calibration/Generated000000.T000.D000.P000.H000.LA.TIF']));
C = imgaussfilt(C,1.5);
c0=[388.5 594.5]; th=(0:1:359)'*pi/180; r=(150:0.25:230)';
pk=nan(size(th)); outer=nan(size(th)); inner=nan(size(th));
for i=1:numel(th)
  xs=c0(1)+r*cos(th(i)); ys=c0(2)-r*sin(th(i)); if any(xs<1|xs>1392|ys<1|ys>1040), continue; end
  p=interp2(C,xs,ys); [m,j]=max(p); bg=median(p(r>215));
  if m-bg<40, continue; end
  pk(i)=r(j); h=bg+(m-bg)/2;
  jo=find(p(j:end)<h,1,'first'); if ~isempty(jo), outer(i)=r(j+jo-1); end
  ji=find(p(1:j)<h,1,'last'); if ~isempty(ji), inner(i)=r(ji); end
end
names={'ring inner half-max','ring peak','ring outer half-max'}; V={inner,pk,outer};
for q=1:3
  rr=V{q}; ok=~isnan(rr);
  for it=1:5
    X=c0(1)+rr(ok).*cos(th(ok)); Y=c0(2)-rr(ok).*sin(th(ok));
    p=[X Y ones(size(X))]\(-(X.^2+Y.^2)); xc=-p(1)/2; yc=-p(2)/2; R=sqrt(xc^2+yc^2-p(3));
    e=hypot(X-xc,Y-yc)-R; s=1.4826*median(abs(e)); idx=find(ok); ok(idx(abs(e)>2.5*s))=false;
  end
  fprintf('%-20s: xc=%.1f row=%.1f D=%.1f px, N=%d, rms=%.2f -> %.5f mm/px\n',names{q},xc,yc,2*R,sum(ok),rms(e),32.15/(2*R));
end
