import type {SidebarsConfig} from '@docusaurus/plugin-content-docs';

const sidebars: SidebarsConfig = {
  docsSidebar: [
    {
      type: 'doc',
      id: 'overview',
      label: 'Overview',
    },
    {
      type: 'category',
      label: 'Get Started',
      items: [
        {
          type: 'doc',
          id: 'usage/quickstart',
          label: 'Quickstart',
        },
        {
          type: 'category',
          label: 'Install',
          collapsed: true,
          items: [
            {
              type: 'doc',
              id: 'installation/docker',
              label: 'Docker',
            },
            {
              type: 'doc',
              id: 'installation/apptainer',
              label: 'Apptainer / Singularity',
            },
            {
              type: 'doc',
              id: 'installation/google-cloud-docker',
              label: 'Google Cloud',
            },
            {
              type: 'doc',
              id: 'installation/non-containerized',
              label: 'From Source',
            },
          ],
        },
        {
          type: 'doc',
          id: 'usage/end-to-end-example-wes',
          label: 'WES/WGS Example',
        },
        {
          type: 'doc',
          id: 'usage/end-to-end-example-mit',
          label: 'mtDNA Example',
        },
        {
          type: 'doc',
          id: 'usage/usage',
          label: 'CLI Usage',
        },
      ],
    },
    {
      type: 'category',
      label: 'Workflows',
      items: [
        {
          type: 'doc',
          id: 'pipelines/overview',
          label: 'Supported Workflows',
        },
        {
          type: 'doc',
          id: 'backends/nf-core',
          label: 'External nf-core',
        },
        {
          type: 'category',
          label: 'Pipeline Guides',
          collapsed: true,
          items: [
            {
              type: 'doc',
              id: 'pipelines/wes-wgs-single',
              label: 'WES/WGS Single-Sample',
            },
            {
              type: 'doc',
              id: 'pipelines/wes-wgs-cohort',
              label: 'WES/WGS Cohort',
            },
            {
              type: 'doc',
              id: 'pipelines/mtdna',
              label: 'mtDNA',
            },
          ],
        },
      ],
    },
    {
      type: 'category',
      label: 'Reproducibility',
      items: [
        {
          type: 'doc',
          id: 'validation/overview',
          label: 'Overview',
        },
        {
          type: 'doc',
          id: 'validation/run-comparison',
          label: 'Run Comparison',
        },
        {
          type: 'category',
          label: 'Validation Evidence',
          collapsed: true,
          items: [
            {
              type: 'doc',
              id: 'validation/integration-tests',
              label: 'Integration Tests',
            },
            {
              type: 'doc',
              id: 'validation/cross-environment',
              label: 'Cross-Environment',
            },
            {
              type: 'doc',
              id: 'validation/giab',
              label: 'GIAB Benchmarking',
            },
            {
              type: 'doc',
              id: 'usage/resource-validation',
              label: 'Resource Validation',
            },
          ],
        },
      ],
    },
    {
      type: 'category',
      label: 'Reference',
      collapsed: true,
      items: [
        {
          type: 'doc',
          id: 'help/configuration-reference',
          label: 'Configuration',
        },
        {
          type: 'doc',
          id: 'help/outputs',
          label: 'Outputs',
        },
        {
          type: 'doc',
          id: 'help/naming-conventions',
          label: 'Naming Conventions',
        },
        {
          type: 'doc',
          id: 'help/performance',
          label: 'Performance',
        },
        {
          type: 'doc',
          id: 'help/troubleshooting',
          label: 'Troubleshooting',
        },
        {
          type: 'doc',
          id: 'help/faq',
          label: 'FAQ',
        },
      ],
    },
    {
      type: 'category',
      label: 'Developers',
      collapsed: true,
      items: [
        {
          type: 'doc',
          id: 'technical-details/architecture',
          label: 'Architecture',
        },
        {
          type: 'doc',
          id: 'technical-details/extending-cbicall',
          label: 'Developer Overview',
        },
        {
          type: 'doc',
          id: 'technical-details/adding-a-pipeline',
          label: 'Adding a Pipeline',
        },
        {
          type: 'doc',
          id: 'technical-details/adding-resources',
          label: 'Adding Resources',
        },
        {
          type: 'doc',
          id: 'technical-details/resource-bundle-v1',
          label: 'Resource Bundle v1',
        },
      ],
    },
    {
      type: 'category',
      label: 'About',
      collapsed: true,
      items: [
        {
          type: 'doc',
          id: 'about/about',
          label: 'About',
        },
        {
          type: 'doc',
          id: 'about/citation',
          label: 'Citation',
        },
        {
          type: 'doc',
          id: 'about/disclaimer',
          label: 'Disclaimer',
        },
      ],
    },
  ],
};

export default sidebars;
