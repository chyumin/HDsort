classdef MessageHandler < handle
    properties (SetAccess=protected)
        
    end
    
    properties
        P
        doSend
    end
    
    methods
        %%% ----------------CONSTRUCTOR------------------------------------
        function self = MessageHandler(varargin)
            gc = global_user_config();
            
            P = struct;
            P.doSend = true;
            
            P.receiver = gc.email_address;
            P.sender = gc.email_address;
            
            P.password = gc.password;
            P.username = gc.user_name;
            P.host = 'mail.ethz.ch';
            P.port = '587';
            P = mysort.hdsort.util.parseInputs(P, varargin, 'error');
            self.P = P;
            
            self.doSend = P.doSend;
            
            setpref('Internet','SMTP_Server', self.P.host);
            setpref('Internet','E_mail', self.P.sender);
            setpref('Internet','SMTP_Username', self.P.username);
            setpref('Internet','SMTP_Password', self.P.password);
            
            props = java.lang.System.getProperties;
            props.setProperty( 'mail.smtp.starttls.enable', 'true' );
            props.setProperty( 'mail.smtp.user', self. P.sender);
            props.setProperty( 'mail.smtp.host', self.P.host );
            props.setProperty( 'mail.smtp.port', self.P.port );
            props.setProperty( 'mail.smtp.debug', 'true' );
            props.setProperty( 'mail.smtp.auth', 'true' );
            props.setProperty( 'mail.smtp.socketFactory.port', self.P.port );
            props.setProperty( 'mail.smtp.socketFactory.fallback', 'false' );
            props.remove( 'mail.smtp.socketFactory.class' );
        end
        
        function sendProgressMessage(self, subject, message)
            if ~self.doSend return; end
            
            assert(~isempty(subject), 'A message title is necessary!')
            if nargin < 3
                message = subject;
            end
            sendmail(self.P.receiver, subject, message)
        end
        
    end
end