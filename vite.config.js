import { defineConfig } from 'vite';
import cesium from 'vite-plugin-cesium';

export default defineConfig({
  // Serve from web/ as root so index.html is at /
  root: 'web',
  // Output build to dist/ at project root
  build: {
    outDir: '../dist',
    emptyOutDir: true,
  },
  plugins: [
    // vite-plugin-cesium handles copying Cesium static assets and setting
    // CESIUM_BASE_URL so the Cesium viewer finds its Workers/Assets/Widgets.
    cesium(),
  ],
  // Allow Vite to resolve the WASM package from propagator/pkg/
  resolve: {
    alias: {
      // propagator/pkg is referenced as '../propagator/pkg/...' in main.js
    },
  },
  // Ensure .wasm files are served correctly
  optimizeDeps: {
    exclude: ['propagator'],
  },
  server: {
    port: 3000,
    open: true,
  },
});
