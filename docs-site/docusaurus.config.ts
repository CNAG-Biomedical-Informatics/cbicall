import {themes as prismThemes} from 'prism-react-renderer';
import type {Config} from '@docusaurus/types';
import type * as Preset from '@docusaurus/preset-classic';

const config: Config = {
  title: 'CBIcall Docs',
  tagline: 'CNAG Biomedical Informatics framework for variant calling',
  favicon: 'img/cbicall-logo.png',
  url: 'https://cnag-biomedical-informatics.github.io',
  baseUrl: '/cbicall/',
  organizationName: 'CNAG-Biomedical-Informatics',
  projectName: 'cbicall',
  onBrokenLinks: 'warn',
  markdown: {
    hooks: {
      onBrokenMarkdownLinks: 'warn',
    },
  },
  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },
  presets: [
    [
      'classic',
      {
        docs: {
          sidebarPath: './sidebars.ts',
          routeBasePath: 'docs',
        },
        blog: false,
        theme: {
          customCss: './src/css/custom.css',
        },
      } satisfies Preset.Options,
    ],
  ],
  themes: [
    [
      '@easyops-cn/docusaurus-search-local',
      {
        hashed: true,
        language: ['en'],
        indexDocs: true,
        indexBlog: false,
        docsRouteBasePath: '/docs',
      },
    ],
  ],
  themeConfig: {
    image: 'img/cbicall-logo.png',
    colorMode: {
      respectPrefersColorScheme: true,
    },
    navbar: {
      title: 'CBIcall',
      logo: {
        alt: 'CBIcall',
        src: 'img/cbicall-logo-small-transparent.png',
      },
      items: [
        {
          type: 'docSidebar',
          sidebarId: 'docsSidebar',
          position: 'left',
          label: 'Docs',
        },
        {
          to: '/docs/usage/quickstart',
          label: 'Quick Start',
          position: 'left',
        },
        {
          type: 'dropdown',
          label: 'Reference',
          position: 'left',
          items: [
            {
              to: '/docs/help/configuration-reference',
              label: 'Configuration Reference',
            },
            {
              to: '/docs/help/outputs',
              label: 'Outputs',
            },
            {
              to: '/docs/help/troubleshooting',
              label: 'Troubleshooting',
            },
            {
              to: '/docs/help/faq',
              label: 'FAQ',
            },
          ],
        },
        {
          href: 'https://github.com/CNAG-Biomedical-Informatics/cbicall',
          label: 'GitHub',
          position: 'right',
        },
      ],
    },
    footer: {
      style: 'dark',
      links: [
        {
          title: 'Docs',
          items: [
            {
              label: 'Overview',
              to: '/docs/overview',
            },
            {
              label: 'Quick Start',
              to: '/docs/usage/quickstart',
            },
            {
              label: 'Troubleshooting',
              to: '/docs/help/troubleshooting',
            },
          ],
        },
        {
          title: 'Project',
          items: [
            {
              label: 'Repository',
              href: 'https://github.com/CNAG-Biomedical-Informatics/cbicall',
            },
            {
              label: 'CNAG',
              href: 'https://www.cnag.eu',
            },
          ],
        },
      ],
      copyright: 'Copyright © 2025-2026 Manuel Rueda, CNAG.',
    },
    prism: {
      theme: prismThemes.github,
      darkTheme: prismThemes.dracula,
    },
  } satisfies Preset.ThemeConfig,
};

export default config;
