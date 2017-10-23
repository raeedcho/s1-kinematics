function unwrapped = remove_wrapping(wrapped)
% Function to remove wrapping from subtracting two angles

wrapped(wrapped>pi) = wrapped(wrapped>pi)-2*pi;
wrapped(wrapped<-pi) = wrapped(wrapped<-pi)+2*pi;
unwrapped = wrapped;

end
