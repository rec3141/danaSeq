import { defineConfig } from 'vite';
import { svelte } from '@sveltejs/vite-plugin-svelte';
import fs from 'fs';
import path from 'path';

// VIZ_DATA_DIR env var overrides public/data for dev/preview serving.
// Usage: VIZ_DATA_DIR=/tmp/viz/my_run npm run dev
const dataDir = process.env.VIZ_DATA_DIR;

export default defineConfig({
  plugins: [
    svelte(),
    // Serve /data/* from VIZ_DATA_DIR when set (avoids polluting source tree)
    ...(dataDir ? [{
      name: 'viz-data-dir',
      configureServer(server) {
        server.middlewares.use('/data', (req, res, next) => {
          const filePath = path.join(dataDir, req.url.split('?')[0]);
          if (fs.existsSync(filePath)) {
            const ext = path.extname(filePath);
            const types = { '.json': 'application/json', '.gz': 'application/json' };
            res.setHeader('Content-Type', types[ext] || 'application/octet-stream');
            if (ext === '.gz') res.setHeader('Content-Encoding', 'gzip');
            fs.createReadStream(filePath).pipe(res);
          } else {
            next();
          }
        });
      },
      configurePreviewServer(server) {
        server.middlewares.use('/data', (req, res, next) => {
          const filePath = path.join(dataDir, req.url.split('?')[0]);
          if (fs.existsSync(filePath)) {
            const ext = path.extname(filePath);
            const types = { '.json': 'application/json', '.gz': 'application/json' };
            res.setHeader('Content-Type', types[ext] || 'application/octet-stream');
            if (ext === '.gz') res.setHeader('Content-Encoding', 'gzip');
            fs.createReadStream(filePath).pipe(res);
          } else {
            next();
          }
        });
      },
    }] : []),
  ],
  server: {
    port: 5173,
    host: '0.0.0.0',
  },
  preview: {
    port: 5174,
    host: '0.0.0.0',
  },
});
