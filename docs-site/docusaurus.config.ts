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
  onBrokenMarkdownLinks: 'warn',
  markdown: {
    mermaid: true,
  },
  themes: ['@docusaurus/theme-mermaid'],
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
  themeConfig: {
    image: 'img/cbicall-logo.png',
    colorMode: {
      respectPrefersColorScheme: true,
    },
    navbar: {
      title: 'CBIcall',
      logo: {
        alt: 'CBIcall',
        src: 'img/cbicall-logo.png',
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
          to: '/docs/usage/choose-your-path',
          label: 'Choose Path',
          position: 'left',
        },
        {
          to: '/docs/help/faq',
          label: 'FAQ',
          position: 'left',
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
