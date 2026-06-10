"""Servidor de desarrollo sin caché para StructWeb3D.
Uso: python serve.py
"""
import http.server, socketserver, os

PORT = 8765

class NoCacheHandler(http.server.SimpleHTTPRequestHandler):
    def guess_type(self, path):
        ctype = super().guess_type(path)
        if isinstance(ctype, str):
            if ctype in ('text/javascript', 'application/javascript', 'text/css', 'text/html'):
                ctype = ctype + '; charset=UTF-8'
        return ctype

    def end_headers(self):
        self.send_header('Cache-Control', 'no-store, no-cache, must-revalidate')
        self.send_header('Pragma', 'no-cache')
        super().end_headers()

    def log_message(self, fmt, *args):
        pass  # silenciar logs de acceso

os.chdir(os.path.dirname(os.path.abspath(__file__)))
with socketserver.TCPServer(('', PORT), NoCacheHandler) as httpd:
    print(f'StructWeb3D → http://localhost:{PORT}')
    httpd.serve_forever()
