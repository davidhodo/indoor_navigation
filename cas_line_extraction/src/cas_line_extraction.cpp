#include "cas_line_extraction/cas_line_extraction.h"
using namespace cas_line_extraction;


LineExtractor::LineExtractor() {

}

LineExtractor::~LineExtractor() {

}

void LineExtractor::ExtractLines(VectorXd scan) {

	unsigned int num_of_scan_points;

//  stdrho = ones(length(phi),1)*stdrho;
//  scan = [phi, rho, stdrho];
//  l = size(scan,1);
//
//
//  % ---------------------------------------------------------- %
//  %%%%%%%%%%% Slide Window And Fit Lines --> alpharC %%%%%%%%%%%
//  % ---------------------------------------------------------- %
//  halfWin = (winsize-1)/2;
//  ibeg = 1 + ~cyc*halfWin;
//  iend = l - ~cyc*halfWin;
//  windowvec = ibeg:iend;
//  for imid = windowvec,
//    for j = imid-halfWin:imid+halfWin,
//      k = mod(j-1,l) + 1;
//      windowpts(j-imid+halfWin+1,:) = scan(k,:);
//    end;
//    [p,C] = fitlinepolar(windowpts);
//    ifillin = mod(imid-1,l) + 1;
//    alpharC(1,ifillin) = p(1,1);			 % alpha
//    alpharC(2,ifillin) = p(2,1);		   % r
//    alpharC(3,ifillin) = C(1,1);		   % sigma aa
//    alpharC(4,ifillin) = C(2,2);			 % sigma rr
//    alpharC(5,ifillin) = C(1,2);	     % sigma ar
//  end;
//
//  % --------------------------------------------------------------------------- %
//  %%%%%%%% Calculate Model Fidelity (=Compactness) Measure -> compactness %%%%%%%
//  % --------------------------------------------------------------------------- %
//  ibeg = 1 + ~cyc*(1+halfWin);   % taking always 3 adjacent points in model space
//  iend = l - ~cyc*(1+halfWin);	 % keep track of shrinking valid range at array ends if cyc=1
//  compactvec  = ibeg:iend;
//  windowlines = zeros(5,3);
//  compactness = zeros(1,l);
//  for imid = compactvec,
//    for j = imid-1:imid+1,
//      k = mod(j-1,l) + 1;
//      windowlines(:,j-imid+2) =  alpharC(:,k);
//    end;
//    compactness(mod(imid-1,l)+1) = calccompactness(windowlines);
//  end;
//
//  % ------------------------------------------------------------------- %
//  %%%%%%%%% Apply Model Fidelity Threshold --> rawSegs(1/2/3,:) %%%%%%%%%
//  % ------------------------------------------------------------------- %
//  nSegBlowup  = halfWin;
//  rSegIndices = findregions(compactness(compactvec),threshfidel,'<',cyc);
//  nRawSegs    = size(rSegIndices,1);
//
//  if rSegIndices(1,1) >= 0,		% check case rSegIndices = -1, see help findregions
//    trwseg = cputime;
//    % Test here for wraparound segment kernels
//    if cyc & (nRawSegs > 1),														        % ingle segm. requires no cyclic treatment
//      o1 = (rSegIndices(1,3) < rSegIndices(1,2));							  % does first segment wrap around?
//      o2 = (rSegIndices(nRawSegs,3) < rSegIndices(nRawSegs,2));	% does last segment wrap around?
//      if (o1 & o2) | (rSegIndices(1,2)-rSegIndices(nRawSegs,3) + ~(o1 | o2)*l <= 1),
//        rSegIndices(1,2) = rSegIndices(nRawSegs,2);             % ibeg of last one is ibeg of new cyclic segm.
//        rSegIndices = rSegIndices(1:nRawSegs-1,:);						  % delete last (=first) raw segment
//        nRawSegs = nRawSegs - 1;
//      end;
//    end;
//    % Blow up segments. Note: if cyc=0 add (1+halfWin) to all indices since they are valid...
//    for i = 1:nRawSegs,																	        % ...for the *shrinked* cmpct-vector
//      rSegIndices(i,2) = mod2(rSegIndices(i,2)-nSegBlowup+~cyc*(1+halfWin),l);
//      rSegIndices(i,3) = mod2(rSegIndices(i,3)+nSegBlowup+~cyc*(1+halfWin),l);
//    end;
//    rawSegs = rSegIndices';
//
//    % --------------------------------------------------------------------------------- %
//    %%%%%%%%% Line Fit To Points In Homogenous Regions -> rawSegs(4/5/6/7/8/9,:) %%%%%%%%
//    % --------------------------------------------------------------------------------- %
//    for i = 1:nRawSegs,
//      segPoints = 0;
//      ibeg = rawSegs(2,i);
//      iend = rawSegs(3,i);
//      if iend < ibeg,                   % discontinuity lies within the segment
//        segPoints = [scan(ibeg:l,:); scan(1:iend, :)];
//      else															% normal case: discontinuity lies not within the segment
//        segPoints = scan(ibeg:iend, :);
//      end;
//      [p,C] = fitlinepolar(segPoints);
//      rawSegs(4,i) = length(segPoints);	% number of raw segment points
//      rawSegs(5,i) = p(1,1);						% alpha
//      rawSegs(6,i) = p(2,1);						% r
//      rawSegs(7,i) = C(1,1);						% sigma_aa
//      rawSegs(8,i) = C(2,2);						% sigma_rr
//      rawSegs(9,i) = C(2,1);						% sigma_ar
//      xybeg = [scan(ibeg,2)*cos(scan(ibeg,1)), scan(ibeg,2)*sin(scan(ibeg,1))];
//      xyend = [scan(iend,2)*cos(scan(iend,1)), scan(iend,2)*sin(scan(iend,1))];
//      Pe1 = calcendpoint(xybeg,rawSegs(5:6,i));
//      Pe2 = calcendpoint(xyend,rawSegs(5:6,i));
//      rawSegs(10:13,i) = [Pe1, Pe2]';	  % Cartesian endpoint coordinates of raw segment
//    end;
//
//    % --------------------------------------------------------------------------- %
//    %%%%%%%%  Agglomerative Hierarchical Clustering --> D,M,rawSegs,nLines %%%%%%%%
//    % --------------------------------------------------------------------------- %
//    nLines = nRawSegs;
//    if nRawSegs > 1,
//      % Fill in distance matrix and set diagonal elements to zero
//      for i=1:nRawSegs-1,
//        for j=i+1:nRawSegs,
//          D(i,j) = mahalanobisar(rawSegs(5:9,i),rawSegs(5:9,j));
//        end;
//      end;
//      for i=1:nRawSegs, D(i,i) = 0; end;
//      M = zeros(nRawSegs);		  % M holds the markers, whether a min in D is valid or not
//
//      ahcmode  = 2;             % mode 2: coll,neigh,over; mode 1: coll,neigh; mode 0: coll
//      terminate = 0;
//      while ~terminate,
//        % Find minimum
//        dmin = Inf;
//        for i = 1:nRawSegs-1,   % finds always first (lowest index) instance of segment in question
//          for j = i+1:nRawSegs,
//            if ~((D(i,j) < 0) | ((ahcmode==1)&(M(i,j)==1)) | ((ahcmode==2)&(M(i,j)==2))),
//              if D(i,j) < dmin,
//                dmin = D(i,j);
//                minrow = i;     % row in D which holds the current minimal distance
//                mincol = j;     % colon in D which holds the current minimal distance
//              end;
//            end;
//          end;
//        end;
//
//        if dmin < fuselevel,
//          if isneighbour(rawSegs, minrow, mincol, ahcmode-1, cyc) | (ahcmode==0),
//            % Fusion of identified segments in 'rawSegs'
//            nLines = nLines - 1;
//            rejectClustNr = rawSegs(1,mincol);
//            for i = 1:nRawSegs,
//              if rawSegs(1,i) == rejectClustNr,
//                rawSegs(1,i) = rawSegs(1,minrow); % mark all segments by overwriting their old ID
//              end;
//            end;
//            for i = 1:nRawSegs,
//              if rawSegs(1,i) == rawSegs(1,minrow),
//                % First: mark all raw points which belong to the segments in question -> scan(:,5)
//                % Simple, less efficient but general solution.
//                ibeg = rawSegs(2,i);
//                iend = rawSegs(3,i);
//                if iend < ibeg,
//                  scan(ibeg:l,5) = rawSegs(1,minrow);
//                  scan(1:iend,5) = rawSegs(1,minrow);
//                else
//                  scan(ibeg:iend, 5) = rawSegs(1,minrow);
//                end;
//              end;
//            end;
//            % Second: go through the scan and collect them -> segPoints
//            segPoints = [-1 -1 -1];
//            segPointCounter = 0;
//            for i = 1:l,
//              if scan(i,5) == rawSegs(1,minrow),
//                segPointCounter = segPointCounter + 1;
//                segPoints = [segPoints; scan(i,1:3)];		% accumulate points
//              end;
//            end;
//            % Third: fit the line with the combined set of points -> new cluster center
//            [p,C] = fitlinepolar(segPoints(2:segPointCounter+1,:));		% remove -1-init-vector at index 1.
//            rawSegs(5:9,minrow) = [p(1,1) p(2,1) C(1,1) C(2,2) C(1,2)]';
//            rawSegs(5:9,mincol) = [p(1,1) p(2,1) C(1,1) C(2,2) C(1,2)]';
//
//            % Calculate new distance in minrow-th row and column
//            for i = minrow+1:nRawSegs,
//              if D(minrow,i) > 0,
//                D(minrow,i) = mahalanobisar(rawSegs(5:9,minrow),rawSegs(5:9,i));
//                M(minrow,i) = 0;   % reset isneighbourhip, v.4.2
//              end;
//            end;
//            for i = 1:minrow-1,
//              if D(i,minrow) > 0,
//                D(i,minrow) = mahalanobisar(rawSegs(5:9,minrow),rawSegs(5:9,i));
//                M(i,minrow) = 0;   % reset isneighbourhip, v.4.2
//              end;
//            end;
//            % Remove mincol-th row and column by marking them with -1
//            for i = mincol+1:nRawSegs,
//              D(mincol,i) = -1; M(mincol,i) = 0;
//            end;	% mark D and unmark M
//            for i = 1:mincol-1,
//              D(i,mincol) = -1; M(i,mincol) = 0;
//            end;	% mark D and unmark M
//          else
//            M(minrow, mincol) = ahcmode; % mark minrow, mincol
//          end;
//        else	% dmin >= fuselevel
//          if ahcmode == 2,
//            if mode == 1,               % mode is changed only if line mode
//              ahcmode = 1;              % is enabled...
//            else
//              terminate = 1;            % ... or terminate otherwise
//            end;
//          elseif ahcmode == 1,
//            ahcmode = 0;                % mode change 2 -> 1
//          else
//            terminate = 1;              % terminate if there are no collinear segments anymore
//          end;
//        end;
//      end;
//    end; % if nRawSegs > 1
//
//    % ---------------------------------------------------------- %
//    %%%%%%%%%% Generate Data Structures --> segs, lines %%%%%%%%%%
//    % ---------------------------------------------------------- %
//    j = 1;
//    for i = min(rawSegs(1,:)):max(rawSegs(1,:));
//      idxvec = find(rawSegs(1,:)==i);
//      if ~isempty(idxvec),
//        nRawSegs = sum(idxvec~=0);						% segIndices holds:
//        segIndices(j,1:length(idxvec)+2) = [i nRawSegs idxvec]; % 1st column: ID of line in 'rawSegs'
//        j = j + 1;														% 2nd column: number of raw segments
//      end;																		% 3:n column: indices of raw segm. in 'rawSegs'
//    end;
//    li = 0; si = 0;
//    for i = 1:nLines, 												% for all lines
//      lineID = segIndices(i,1); 							% get ID
//      for j = 1:segIndices(i,2);              % for all raw segments constituting the line lineID
//        ibeg = rawSegs(2,segIndices(i,2+j));  % get their begin- and end-indices which,
//        iend = rawSegs(3,segIndices(i,2+j));  % in general, overlap the adjacent segments
//        if iend < ibeg, 										  % discontinuity lies within segment
//          scan(ibeg:l,5) = lineID; 						% mark raw points
//          scan(1:iend,5) = lineID;
//        else  															  % normal case: discontinuity lies not within segment
//          scan(ibeg:iend,5) = lineID;         % mark raw points
//        end;
//      end;
//      nLinePoints = (sum(scan(:,5) == lineID));
//      jSegIndices = findregions(scan(:,5),lineID,'==',cyc);		% find joint segments
//      nJointSegs  = size(jSegIndices, 1);
//
//      lmax = -Inf; imax = 0;
//      for j = 1:nJointSegs,
//        ibeg = jSegIndices(j,2); iend = jSegIndices(j,3);
//        [pBeg(1),pBeg(2)] = pol2cart(scan(ibeg,1),scan(ibeg,2));
//        [pEnd(1),pEnd(2)] = pol2cart(scan(iend,1),scan(iend,2));
//        Pe1 = calcendpoint(pBeg,rawSegs(5:6,segIndices(i,3))); % Cartesian projection
//        Pe2 = calcendpoint(pEnd,rawSegs(5:6,segIndices(i,3)));
//        len = sqrt(sum((Pe1-Pe2).^2));
//        if len > lmax,
//          lmax = len; imax = j;
//        end;
//      end;
//
//      if lmax >= minlength,
//        li = li + 1;
//        lines(1,li) = lineID;                     % line ID
//        lines(2,li) = nJointSegs;                 % number of contributing joint segments
//        lines(3,li) = segIndices(i,2);				    % number of contributing raw segments
//        lines(4,li) = nLinePoints;                % number of contributing raw data points
//        lines(5:6,li) = rawSegs(5:6,segIndices(i,3));           % alpha,r
//        lines(7,li) = rawSegs(7,segIndices(i,3)) + saacompens;  % saa
//        lines(8,li) = rawSegs(8,segIndices(i,3)) + srrcompens;  % srr
//        lines(9,li) = rawSegs(9,segIndices(i,3));               % sar
//        lines(10,li) = sum(rawSegs(7:8,segIndices(i,3))); % trace of covariance matrix
//        for j = 1:nJointSegs,
//          si = si + 1;
//          lines(10+j,li) = si; 				  	        % IDs of contributing joint segments
//          segs(1,si) = si; 					              % ID of joint segment
//          segs(2,si) = lineID; 					          % ID of line the segment belongs to
//          segs(3,si) = segIndices(i,2);	          % number of contributing raw segments
//          nPoints = mod2(jSegIndices(j,3)-jSegIndices(j,2)+1, l);
//          segs(4,si) = nPoints;					          % number of contributing raw data points
//          segs(5:10, si) = lines(5:10,li);				% copy alpha,r and cov-matrix
//          segs(11:12,si) = jSegIndices(j,2:3)';   % begin, end index of joint segment
//          ibeg = jSegIndices(j,2);
//          iend = jSegIndices(j,3);
//          [pBeg(1),pBeg(2)] = pol2cart(scan(ibeg,1), scan(ibeg,2));
//          [pEnd(1),pEnd(2)] = pol2cart(scan(iend,1), scan(iend,2));
//          Pe1 = calcendpoint(pBeg,lines(5:6,li)); % Cartesian projection
//          Pe2 = calcendpoint(pEnd,lines(5:6,li));
//          segs(13:16,si) = [Pe1, Pe2]';	          % Cartesian coordinates of endpoints
//        end;
//      end;
//    end;
//    nJointSegs = si;
//    nLines = li;
//
//    % Resort 'segs' such that it is sorted acc. to real succession of segments in scan
//    if nLines > 0
//      [dummy,jSegsisort] = sortrows(segs',11);
//      for i = 1:nJointSegs,
//        tmp(:,i) = segs(:,jSegsisort(i));
//      end;
//      segs = tmp;
//
//      % ---------------------------------------------------------- %
//      %%%%%%%%%%%  Generate Final Data Structures --> L  %%%%%%%%%%%
//      % ---------------------------------------------------------- %
//
//      % Put everything into the new map object L
//      created = 0; iline = 1;
//      for i = 1:nLines;
//        % create map if not yet done
//        if ~created,
//          L = map('local map',0);
//          created = 1;
//        end;
//        % get alpha-r values
//        id    = lines(1,i);
//        alpha = lines(5,i);
//        r     = lines(6,i);
//        saa   = lines(7,i);
//        srr   = lines(8,i);
//        sar   = lines(9,i);
//        % create line feature
//        l = arlinefeature(id,[alpha;r],[saa sar; sar srr]);
//        % add it to local map
//        L = addentity(L,l);
//      end;
//      if created,
//        X = get(L,'x');
//        nL = length(X);
//      else
//        L = []; nL = 0;
//      end;
//
//      if nL == 0, str1 = 'no'; str2 = 'lines';
//      elseif nL == 1, str1 = '1'; str2 = 'line';
//      else str1 = int2str(nL); str2 = 'lines';
//      end; %disp([' ',str1,' ',str2,' extracted']);
//
//      % Assign output
//      if nargout == 1;
//        varargout{1} = L;
//      else
//        varargout{1} = L;
//        varargout{2} = segs;
//        varargout{3} = lines;
//      end;




}
