
function logToFile(filename, str)
  hdsort.util.appendToFile(filename, [datestr(now) '  ' str '\n']);
end